"""
make_viewer_gmx.py — Standalone HTML trajectory viewer for GROMACS MD output.
Works directly from file:// — no server needed.

Usage:
    python make_viewer_gmx.py [--outdir results/run]
    python make_viewer_gmx.py --traj md.xtc --top md.tpr --outdir .
"""

import argparse
import os
import sys
from pathlib import Path

# Inline JS libraries (no CDN needed — works offline / file://)
_SCRIPT_DIR = Path(__file__).parent
_JQUERY_JS = (_SCRIPT_DIR / "jquery.min.js").read_text()
_3DMOL_JS  = (_SCRIPT_DIR / "3Dmol-min.js").read_text()

import mdtraj as md
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def main():
    parser = argparse.ArgumentParser(description="Generate HTML viewer from GROMACS trajectory")
    parser.add_argument("--outdir", default=".", help="Directory with md.xtc / md.tpr (and output)")
    parser.add_argument("--traj",   default=None)
    parser.add_argument("--top",    default=None)
    parser.add_argument("--ns",     type=float, default=None, help=argparse.SUPPRESS)  # deprecated, auto-detected
    parser.add_argument("--stride", type=int,   default=10,
                        help="Load every Nth frame (default 10). Increase for large trajectories.")
    args = parser.parse_args()

    outdir = Path(args.outdir)
    traj_p = Path(args.traj) if args.traj else outdir / "md.xtc"
    top_p  = Path(args.top)  if args.top  else outdir / "npt.gro"

    for p in (traj_p, top_p):
        if not p.exists():
            sys.exit(f"File not found: {p}")

    # ── Load & strip solvent ────────────────────────────────────────────────
    print(f"Loading trajectory: {traj_p}  (stride={args.stride})")
    traj = md.load(str(traj_p), top=str(top_p), stride=args.stride)
    print(f"  Raw: {traj.n_frames} frames, {traj.n_atoms} atoms")

    traj.image_molecules(inplace=True)

    solute_idx = traj.topology.select(
        "not (water or resname SOL or resname HOH or resname NA or resname CL)"
    )
    traj = traj.atom_slice(solute_idx)
    traj.superpose(traj, 0)
    traj_full = traj  # keep full trajectory for plots

    # Downsample to ≤100 frames for browser performance
    max_frames = 100
    if traj.n_frames > max_frames:
        stride = traj.n_frames // max_frames
        traj = traj[::stride]
    print(f"  Stripped: {traj.n_atoms} solute atoms, {traj.n_frames} frames (full: {traj_full.n_frames})")

    # ── Detect ligand residue name ──────────────────────────────────────────
    lig_resname = None
    for res in traj.topology.residues:
        if not res.is_protein and res.name not in ("HOH", "SOL", "WAT", "NA", "CL"):
            lig_resname = res.name
            break
    if lig_resname is None:
        lig_resname = "MOL"
    print(f"  Ligand residue name: {lig_resname}")

    # ── Names for title ─────────────────────────────────────────────────────
    folder = outdir.resolve().name          # e.g. 4aqp_chainA_receptor_top_pose_1_GMX
    parts  = folder.replace("_GMX", "").split("_")
    # heuristic: split at first part that looks like a compound name
    pdb_name      = folder.split("_receptor")[0] if "_receptor" in folder else parts[0]
    compound_name = folder.split("_GMX")[0].split("_")[-2] + "_" + folder.split("_GMX")[0].split("_")[-1] \
                    if "_" in folder else "ligand"
    title = f"{pdb_name}  ·  {compound_name}"

    # ── Compute inline metrics ──────────────────────────────────────────────
    ca_idx     = traj.topology.select("protein and name CA")
    lig_idx    = traj.topology.select(f"resname {lig_resname} and not element H")

    ref = traj[0]
    traj.superpose(ref, atom_indices=ca_idx)
    protein_rmsd = (md.rmsd(traj, ref, 0, atom_indices=ca_idx, precentered=True) * 10).tolist()
    ligand_rmsd  = (md.rmsd(traj, ref, 0, atom_indices=lig_idx, precentered=True) * 10).tolist() if len(lig_idx) else []

    all_lig  = set(traj.topology.select(f"resname {lig_resname}").tolist())
    all_prot = set(traj.topology.select("protein").tolist())

    polar_elements = {"N", "O", "S"}
    lig_polar  = np.array([a.index for a in traj.topology.atoms
                           if a.residue.name == lig_resname and a.element.symbol in polar_elements])
    prot_polar = np.array([a.index for a in traj.topology.atoms
                           if a.residue.is_protein and a.element.symbol in polar_elements])

    hbond_counts   = []
    contact_counts = []
    for i in range(traj.n_frames):
        hb = md.baker_hubbard(traj[i], periodic=False)
        hbond_counts.append(sum(
            1 for d, h, a in hb
            if (int(d) in all_prot and int(a) in all_lig) or
               (int(d) in all_lig  and int(a) in all_prot)
        ))
        if len(lig_polar) > 0 and len(prot_polar) > 0:
            pairs = np.array([[lp, pp] for lp in lig_polar for pp in prot_polar])
            dists = md.compute_distances(traj[i], pairs)[0]
            contact_counts.append(int((dists < 0.35).sum()))
        else:
            contact_counts.append(0)

    # ── Generate 4 analysis plots ───────────────────────────────────────────
    print("Generating plots...")
    t_ns = traj.time / 1000.0

    # Plot A — H-bonds + polar contacts
    hb_arr = np.array(hbond_counts)
    ct_arr = np.array(contact_counts)
    fig, ax = plt.subplots(figsize=(10, 4))
    ax.plot(t_ns, ct_arr, lw=0.8, color="steelblue", alpha=0.5,
            label=f"Polar contacts ≤3.5 Å  (mean {ct_arr.mean():.2f})")
    ax.plot(t_ns, hb_arr, lw=0.8, color="darkorange", alpha=0.8,
            label=f"H-bonds (strict)  (mean {hb_arr.mean():.2f})")
    ax.axhline(ct_arr.mean(), color="steelblue", lw=1.2, ls="--", alpha=0.7)
    ax.axhline(hb_arr.mean(), color="crimson",   lw=1.5, ls="--",
               label=f"Mean H-bonds = {hb_arr.mean():.2f}")
    ax.set_xlabel("Time (ns)"); ax.set_ylabel("# interactions")
    ax.set_title("A — Protein–Ligand Hydrogen Bonds & Polar Contacts"); ax.legend(fontsize=8)
    fig.tight_layout(); fig.savefig(outdir / "plot_A_hbonds.png", dpi=150); plt.close(fig)
    print(f"  → plot_A_hbonds.png")

    # Plot B — RMSD (full traj for smoother curves)
    ca_all  = traj_full.topology.select("protein and name CA")
    lig_all = traj_full.topology.select(f"resname {lig_resname} and not element H")
    ref_full = traj_full[0]
    traj_full.superpose(ref_full, atom_indices=ca_all)
    rmsd_p  = md.rmsd(traj_full, ref_full, 0, atom_indices=ca_all, precentered=True) * 10
    rmsd_l  = md.rmsd(traj_full, ref_full, 0, atom_indices=lig_all, precentered=True) * 10 if len(lig_all) else None
    t_full  = traj_full.time / 1000.0
    fig, ax = plt.subplots(figsize=(10, 4))
    ax.plot(t_full, rmsd_p, lw=0.8, label="Protein Cα", color="steelblue")
    if rmsd_l is not None:
        ax.plot(t_full, rmsd_l, lw=0.8, label=f"Ligand ({lig_resname})", color="darkorange")
    ax.set_xlabel("Time (ns)"); ax.set_ylabel("RMSD (Å)"); ax.set_title("B — RMSD from Starting Structure"); ax.legend()
    fig.tight_layout(); fig.savefig(outdir / "plot_B_rmsd.png", dpi=150); plt.close(fig)
    print(f"  → plot_B_rmsd.png")

    # Plot C — RMSF
    rmsf   = md.rmsf(traj_full, traj_full, 0, atom_indices=ca_all) * 10
    res_ids = [traj_full.topology.atom(i).residue.resSeq for i in ca_all]
    fig, ax = plt.subplots(figsize=(12, 4))
    ax.plot(res_ids, rmsf, lw=0.8, color="mediumseagreen")
    ax.fill_between(res_ids, rmsf, alpha=0.3, color="mediumseagreen")
    ax.set_xlabel("Residue number"); ax.set_ylabel("RMSF (Å)"); ax.set_title("C — Per-Residue Cα RMSF")
    fig.tight_layout(); fig.savefig(outdir / "plot_C_rmsf.png", dpi=150); plt.close(fig)
    print(f"  → plot_C_rmsf.png")

    # Plot D — Rg
    prot_idx = traj_full.topology.select("protein")
    rg = md.compute_rg(traj_full.atom_slice(prot_idx)) * 10
    fig, ax = plt.subplots(figsize=(10, 4))
    ax.plot(t_full, rg, lw=0.8, color="mediumpurple")
    ax.axhline(rg.mean(), color="crimson", lw=1.5, ls="--", label=f"Mean = {rg.mean():.2f} Å")
    ax.set_xlabel("Time (ns)"); ax.set_ylabel("Rg (Å)"); ax.set_title("D — Protein Radius of Gyration"); ax.legend()
    fig.tight_layout(); fig.savefig(outdir / "plot_D_rg.png", dpi=150); plt.close(fig)
    print(f"  → plot_D_rg.png")

    # ── Save stripped PDB for inline embedding ──────────────────────────────
    pdb_path = outdir / "trajectory_solute_gmx.pdb"
    traj.save_pdb(str(pdb_path))
    pdb_data = pdb_path.read_text()
    pdb_js   = pdb_data.replace('\\', '\\\\').replace('`', '\\`').replace('$', '\\$')

    n_frames = traj.n_frames
    total_ns = round(traj_full.time[-1] / 1000, 3)

    # ── Build HTML (identical layout to OpenMM viewer) ──────────────────────
    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>MD Viewer (GROMACS) — {title}</title>
<script>{_JQUERY_JS}</script>
<script>{_3DMOL_JS}</script>
<style>
  *, *::before, *::after {{ box-sizing: border-box; margin: 0; padding: 0; }}
  :root {{
    --bg:      #0d0f1a;
    --panel:   #13152a;
    --border:  #1e2240;
    --accent:  #4f8ef7;
    --accent2: #a78bfa;
    --text:    #e2e8f0;
    --muted:   #64748b;
    --green:   #34d399;
    --orange:  #fb923c;
  }}
  body {{ background:var(--bg); color:var(--text); font-family:'Inter',system-ui,sans-serif;
          display:flex; flex-direction:column; height:100vh; overflow:hidden; }}
  #header {{ background:var(--panel); border-bottom:1px solid var(--border); padding:10px 20px;
             display:flex; align-items:center; justify-content:space-between; flex-shrink:0; z-index:10; }}
  #header h1 {{ font-size:14px; font-weight:600; letter-spacing:.03em; color:var(--text); }}
  #header h1 span {{ color:var(--accent); }}
  .badge {{ font-size:11px; padding:2px 8px; border-radius:99px; background:var(--border); color:var(--muted); font-weight:500; }}
  .badge.green {{ background:#052e16; color:var(--green); }}
  .badge.gmx   {{ background:#1a0e2a; color:var(--accent2); }}
  #main {{ display:flex; flex:1; overflow:hidden; }}
  #viewer {{ flex:1; position:relative; }}
  #sidebar {{ width:220px; background:var(--panel); border-left:1px solid var(--border);
              display:flex; flex-direction:column; overflow-y:auto; flex-shrink:0; }}
  .panel-section {{ padding:14px 16px; border-bottom:1px solid var(--border); }}
  .panel-section h3 {{ font-size:10px; font-weight:700; letter-spacing:.1em; text-transform:uppercase;
                       color:var(--muted); margin-bottom:10px; }}
  .stat-row {{ display:flex; justify-content:space-between; align-items:baseline; margin-bottom:6px; }}
  .stat-label {{ font-size:11px; color:var(--muted); }}
  .stat-val {{ font-size:13px; font-weight:600; color:var(--text); }}
  .stat-val.accent {{ color:var(--accent); }}
  canvas.spark {{ width:100%; height:48px; border-radius:4px; display:block; margin-top:6px; }}
  .style-btn {{ display:block; width:100%; padding:7px 10px; margin-bottom:6px; background:var(--border);
                border:1px solid transparent; border-radius:6px; color:var(--text); font-size:12px;
                cursor:pointer; text-align:left; transition:border-color .15s,background .15s; }}
  .style-btn:hover  {{ background:#1e2650; }}
  .style-btn.active {{ border-color:var(--accent); color:var(--accent); }}
  #controls {{ background:var(--panel); border-top:1px solid var(--border); padding:10px 20px;
               display:flex; align-items:center; gap:14px; flex-shrink:0; }}
  .ctrl-btn {{ background:var(--border); border:none; color:var(--text); border-radius:6px;
               width:32px; height:32px; cursor:pointer; font-size:14px;
               display:flex; align-items:center; justify-content:center; transition:background .15s; flex-shrink:0; }}
  .ctrl-btn:hover {{ background:var(--accent); }}
  #playBtn {{ width:40px; height:32px; background:var(--accent); color:#fff; font-size:13px; }}
  #playBtn:hover {{ background:#3b82f6; }}
  #timeline-wrap {{ flex:1; display:flex; flex-direction:column; gap:4px; }}
  #frameSlider {{ width:100%; cursor:pointer; accent-color:var(--accent); height:4px; }}
  #timeline-labels {{ display:flex; justify-content:space-between; font-size:10px; color:var(--muted); }}
  .ctrl-group {{ display:flex; align-items:center; gap:8px; flex-shrink:0; }}
  .ctrl-label {{ font-size:11px; color:var(--muted); }}
  #speedSlider {{ width:80px; accent-color:var(--accent2); }}
  #timeDisplay {{ font-size:12px; color:var(--text); font-variant-numeric:tabular-nums; min-width:110px; flex-shrink:0; }}
  #timeDisplay span {{ color:var(--accent); font-weight:600; }}
  #keyhint {{ font-size:10px; color:var(--muted); flex-shrink:0; }}
  kbd {{ background:var(--border); border-radius:3px; padding:1px 5px; font-size:10px; font-family:monospace; }}
</style>
</head>
<body>

<div id="header">
  <h1>MD Trajectory &nbsp;·&nbsp; <span>{pdb_name}</span> + {compound_name}</h1>
  <div style="display:flex;gap:8px;align-items:center">
    <span class="badge gmx">GROMACS</span>
    <span class="badge">{total_ns} ns simulation</span>
    <span class="badge">{n_frames} frames</span>
    <span class="badge green" id="statusBadge">● Stable</span>
  </div>
</div>

<div id="main">
  <div id="viewer"></div>
  <div id="sidebar">

    <div class="panel-section">
      <h3>Current Frame</h3>
      <div class="stat-row"><span class="stat-label">Time</span><span class="stat-val accent" id="sideTime">0.0 ns</span></div>
      <div class="stat-row"><span class="stat-label">Frame</span><span class="stat-val" id="sideFrame">0 / {n_frames-1}</span></div>
      <div class="stat-row"><span class="stat-label">H-bonds</span><span class="stat-val" id="sideHbonds">—</span></div>
      <div class="stat-row"><span class="stat-label">Prot. RMSD</span><span class="stat-val" id="sideRmsd">—</span></div>
      <div class="stat-row"><span class="stat-label">Lig. RMSD</span><span class="stat-val" id="sideLigRmsd">—</span></div>
    </div>

    <div class="panel-section"><h3>Protein RMSD (Å)</h3><canvas class="spark" id="sparkRmsd"></canvas></div>
    <div class="panel-section"><h3>Ligand RMSD (Å)</h3><canvas class="spark" id="sparkLigRmsd"></canvas></div>
    <div class="panel-section"><h3>H-bonds</h3><canvas class="spark" id="sparkHbond"></canvas></div>

    <div class="panel-section">
      <h3>Representation</h3>
      <button class="style-btn active" id="btnCartoon">Cartoon + Ligand</button>
      <button class="style-btn" id="btnSurface">Surface + Ligand</button>
      <button class="style-btn" id="btnWire">Wireframe</button>
    </div>
    <div class="panel-section">
      <h3>View</h3>
      <button class="style-btn" id="btnZoomLig">Focus ligand</button>
      <button class="style-btn" id="btnZoomAll">Zoom to full</button>
    </div>

  </div>
</div>

<div id="controls">
  <button class="ctrl-btn" id="stepBack">&#9664;</button>
  <button class="ctrl-btn" id="playBtn">&#9654;</button>
  <button class="ctrl-btn" id="stepFwd">&#9654;&#9654;</button>
  <div id="timeline-wrap">
    <input type="range" id="frameSlider" min="0" max="{n_frames-1}" value="0">
    <div id="timeline-labels"><span>0 ns</span><span>{total_ns/2:.1f} ns</span><span>{total_ns} ns</span></div>
  </div>
  <div id="timeDisplay">Time: <span id="nsDisplay">0.0</span> ns</div>
  <div class="ctrl-group"><span class="ctrl-label">Speed</span><input type="range" id="speedSlider" min="50" max="1000" value="200"></div>
  <div id="keyhint"><kbd>Space</kbd> play &nbsp; <kbd>←</kbd><kbd>→</kbd> step</div>
</div>

<script>
const PDB           = `{pdb_js}`;
const N_FRAMES      = {n_frames};
const NS_TOTAL      = {total_ns};
const LIG_RESN      = "{lig_resname}";
const RMSD_DATA     = {protein_rmsd};
const LIG_RMSD_DATA = {ligand_rmsd};
const HBOND_DATA    = {hbond_counts};

let viewer = $3Dmol.createViewer("viewer", {{ backgroundColor:"0x0d0f1a", antialias:true }});
viewer.addModelsAsFrames(PDB, "pdb");

function applyCartoon() {{
  viewer.setStyle({{}}, {{ cartoon: {{ color:"spectrum", opacity:0.95 }} }});
  viewer.setStyle({{ resn: LIG_RESN }}, {{
    stick:  {{ colorscheme:"default", radius:0.18 }},
    sphere: {{ colorscheme:"default", radius:0.38 }}
  }});
}}
applyCartoon();
viewer.zoomTo(); viewer.render(); viewer.setFrame(0);

function drawSparkline(canvasId, data, color, fillColor) {{
  const canvas = document.getElementById(canvasId);
  if (!canvas || !data.length) return null;
  const dpr = window.devicePixelRatio || 1;
  canvas.width  = canvas.offsetWidth  * dpr;
  canvas.height = canvas.offsetHeight * dpr;
  const ctx = canvas.getContext('2d');
  ctx.scale(dpr, dpr);
  const W = canvas.offsetWidth, H = canvas.offsetHeight;
  const mn = Math.min(...data), mx = Math.max(...data) || 1;
  const x = i => (i / (data.length - 1)) * W;
  const y = v => H - ((v - mn) / (mx - mn)) * (H - 6) - 3;
  ctx.beginPath(); ctx.moveTo(x(0), y(data[0]));
  data.forEach((v,i) => ctx.lineTo(x(i), y(v)));
  ctx.strokeStyle = color; ctx.lineWidth = 1.5; ctx.stroke();
  ctx.lineTo(x(data.length-1), H); ctx.lineTo(0, H); ctx.closePath();
  ctx.fillStyle = fillColor; ctx.fill();
  return {{ x, y, mn, mx, W, H }};
}}

let rmsdMeta, ligRmsdMeta, hbondMeta;
function initSparklines() {{
  rmsdMeta    = drawSparkline('sparkRmsd',    RMSD_DATA,     '#4f8ef7', 'rgba(79,142,247,0.15)');
  ligRmsdMeta = drawSparkline('sparkLigRmsd', LIG_RMSD_DATA, '#34d399', 'rgba(52,211,153,0.15)');
  hbondMeta   = drawSparkline('sparkHbond',   HBOND_DATA,    '#a78bfa', 'rgba(167,139,250,0.15)');
}}
setTimeout(initSparklines, 200);

function redrawSpark(canvasId, data, color, fill, frame) {{
  const m = drawSparkline(canvasId, data, color, fill);
  if (!m) return;
  const canvas = document.getElementById(canvasId);
  const dpr = window.devicePixelRatio || 1;
  const ctx = canvas.getContext('2d');
  const W = canvas.offsetWidth, H = canvas.offsetHeight;
  const cx = (frame / (data.length - 1)) * W;
  ctx.beginPath(); ctx.moveTo(cx, 0); ctx.lineTo(cx, H);
  ctx.strokeStyle = '#ffffff55'; ctx.lineWidth = 1; ctx.stroke();
}}

let currentFrame = 0, playing = false, speed = 850, playInterval = null, currentStyle = 'cartoon';
let currentSurfaceId = null, surfaceBusy = false;

function nsFromFrame(f) {{ return (f / (N_FRAMES - 1) * NS_TOTAL).toFixed(2); }}

function showFrame(idx) {{
  currentFrame = idx;
  viewer.setFrame(idx);
  if (currentStyle === 'surface' && !surfaceBusy) {{
    surfaceBusy = true;
    const oldId = currentSurfaceId;
    viewer.addSurface('MS', {{ opacity:1.0, colorscheme:'whiteCarbon' }}, {{ not:{{ resn:LIG_RESN }} }}).then(sid => {{
      currentSurfaceId = sid;
      if (oldId !== null) viewer.removeSurface(oldId);
      viewer.render(); surfaceBusy = false;
    }});
  }} else {{ viewer.render(); }}
  document.getElementById('frameSlider').value = idx;
  document.getElementById('nsDisplay').textContent   = nsFromFrame(idx);
  document.getElementById('sideTime').textContent    = nsFromFrame(idx) + ' ns';
  document.getElementById('sideFrame').textContent   = idx + ' / ' + (N_FRAMES-1);
  document.getElementById('sideHbonds').textContent  = HBOND_DATA[idx] ?? '—';
  document.getElementById('sideRmsd').textContent    = (RMSD_DATA[idx] ?? 0).toFixed(2) + ' Å';
  document.getElementById('sideLigRmsd').textContent = LIG_RMSD_DATA.length ? (LIG_RMSD_DATA[idx] ?? 0).toFixed(2) + ' Å' : '—';
  redrawSpark('sparkRmsd',    RMSD_DATA,     '#4f8ef7', 'rgba(79,142,247,0.15)',   idx);
  redrawSpark('sparkLigRmsd', LIG_RMSD_DATA, '#34d399', 'rgba(52,211,153,0.15)',   idx);
  redrawSpark('sparkHbond',   HBOND_DATA,    '#a78bfa', 'rgba(167,139,250,0.15)',  idx);
}}

function togglePlay() {{
  playing = !playing;
  document.getElementById('playBtn').innerHTML = playing ? '&#9646;&#9646;' : '&#9654;';
  if (playing) {{ playInterval = setInterval(() => showFrame((currentFrame + 1) % N_FRAMES), speed); }}
  else {{ clearInterval(playInterval); }}
}}

document.getElementById('playBtn').addEventListener('click', togglePlay);
document.getElementById('stepBack').addEventListener('click', () => {{ if(playing) togglePlay(); showFrame((currentFrame-1+N_FRAMES)%N_FRAMES); }});
document.getElementById('stepFwd').addEventListener('click',  () => {{ if(playing) togglePlay(); showFrame((currentFrame+1)%N_FRAMES); }});
document.getElementById('frameSlider').addEventListener('input', function() {{ showFrame(+this.value); }});
document.getElementById('speedSlider').addEventListener('input', function() {{
  speed = 1050 - +this.value;
  if (playing) {{ clearInterval(playInterval); playInterval = setInterval(() => showFrame((currentFrame+1)%N_FRAMES), speed); }}
}});
document.addEventListener('keydown', e => {{
  if (e.code==='Space')      {{ e.preventDefault(); togglePlay(); }}
  if (e.code==='ArrowRight') {{ if(playing) togglePlay(); showFrame((currentFrame+1)%N_FRAMES); }}
  if (e.code==='ArrowLeft')  {{ if(playing) togglePlay(); showFrame((currentFrame-1+N_FRAMES)%N_FRAMES); }}
}});

function setActive(id) {{ document.querySelectorAll('.style-btn').forEach(b => b.classList.remove('active')); document.getElementById(id).classList.add('active'); }}
document.getElementById('btnCartoon').addEventListener('click', () => {{ currentStyle='cartoon'; viewer.removeAllSurfaces(); applyCartoon(); viewer.render(); setActive('btnCartoon'); }});
document.getElementById('btnSurface').addEventListener('click', () => {{
  currentStyle='surface';
  viewer.setStyle({{}},{{}});
  viewer.setStyle({{ resn:LIG_RESN }}, {{ stick:{{ colorscheme:'default',radius:0.2 }}, sphere:{{ colorscheme:'default',radius:0.4 }} }});
  viewer.addSurface('MS', {{ opacity:1.0, colorscheme:'whiteCarbon' }}, {{ not:{{ resn:LIG_RESN }} }}).then(sid => {{ currentSurfaceId=sid; viewer.render(); }});
  setActive('btnSurface');
}});
document.getElementById('btnWire').addEventListener('click', () => {{ currentStyle='wire'; viewer.removeAllSurfaces(); viewer.setStyle({{}},{{ line:{{}} }}); viewer.setStyle({{ resn:LIG_RESN }},{{ stick:{{ colorscheme:'default',radius:0.15 }} }}); viewer.render(); setActive('btnWire'); }});
document.getElementById('btnZoomLig').addEventListener('click', () => {{ viewer.zoomTo({{ resn:LIG_RESN }}); viewer.render(); }});
document.getElementById('btnZoomAll').addEventListener('click', () => {{ viewer.zoomTo(); viewer.render(); }});
</script>
</body>
</html>"""

    out = outdir / "trajectory_viewer_gmx.html"
    out.write_text(html)
    size_mb = out.stat().st_size / 1e6
    print(f"\nSaved: {out}  ({size_mb:.1f} MB) — open in browser!")


if __name__ == "__main__":
    main()
