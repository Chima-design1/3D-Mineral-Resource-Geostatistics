"""Render a two-dimensional projection of samples and supported 3D grid nodes."""
import csv
import json
import math
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
RAW=ROOT/"data"/"raw"/"drillholes.csv"
GRID=ROOT/"outputs"/"copper_idw_grid_below_proxy_surface.csv"
FIG=ROOT/"outputs"/"portfolio_spatial_figure.svg"
META=ROOT/"reports"/"portfolio_spatial_figure.json"


def read(path):
    with path.open(encoding="utf-8-sig",newline="") as f:return list(csv.DictReader(f))


def color(value,vmax):
    stops=[(0.0,(49,54,149)),(.25,(69,117,180)),(.5,(116,173,209)),(.7,(254,224,139)),(.85,(244,109,67)),(1.0,(165,0,38))]
    t=math.log1p(max(0,value))/math.log1p(vmax) if vmax>0 else 0
    for i in range(1,len(stops)):
        if t<=stops[i][0]:
            t0,c0=stops[i-1];t1,c1=stops[i];f=(t-t0)/(t1-t0)
            rgb=tuple(round(a+(b-a)*f) for a,b in zip(c0,c1));return f"rgb({rgb[0]},{rgb[1]},{rgb[2]})"
    return "rgb(165,0,38)"


def main():
    samples=read(RAW);grid=read(GRID)
    sample=[(float(r["X"]),float(r["Y"]),float(r["Z"]),float(r["Cu ppm"])) for r in samples]
    estimates=[(float(r["X"]),float(r["Y"]),float(r["Z"]),float(r["Cu_IDW_ppm"])) for r in grid]
    all_values=[p[3] for p in sample]+[p[3] for p in estimates];vmax=max(all_values)
    xmin=min(p[0] for p in sample);xmax=max(p[0] for p in sample);ymin=min(p[1] for p in sample);ymax=max(p[1] for p in sample)
    zmin=min(p[2] for p in sample);zmax=max(p[2] for p in sample);cuts=[zmin+(zmax-zmin)/3,zmin+2*(zmax-zmin)/3]
    slices=[([p for p in estimates if p[2]<cuts[0]],f"Lower Z third: Z < {cuts[0]:.0f}"),([p for p in estimates if cuts[0]<=p[2]<cuts[1]],f"Middle Z third: {cuts[0]:.0f} ≤ Z < {cuts[1]:.0f}"),([p for p in estimates if p[2]>=cuts[1]],f"Upper Z third: Z ≥ {cuts[1]:.0f}")]
    panels=[(sample,"Observed sample Cu",2.2)]+[(ps,title,4.0) for ps,title in slices]
    W,H=1200,720;plotw,ploth=470,235;body=[f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}">','<rect width="100%" height="100%" fill="#ffffff"/>','<style>text{font-family:Arial,sans-serif;fill:#263445}.main{font-size:26px;font-weight:700}.sub{font-size:14px;fill:#475569}.pt{stroke:#fff;stroke-width:.35}.axis{stroke:#526174;stroke-width:1}.tick{font-size:11px;fill:#526174}.panel{font-size:16px;font-weight:600}.note{font-size:12px;fill:#475569}</style>']
    body.extend(['<text class="main" x="600" y="38" text-anchor="middle">GeoMet Cu samples and supported IDW grid</text>','<text class="sub" x="600" y="63" text-anchor="middle">XY projections of 3D points • color = log-scaled Cu concentration (ppm) • grid is exploratory, not a resource model</text>'])
    for ix,(points,title,radius) in enumerate(panels):
        col=ix%2;row=ix//2;left=75+col*560;top=95+row*285
        body.append(f'<text class="panel" x="{left}" y="{top-12}">{title} (n={len(points)})</text>')
        body.append(f'<rect x="{left}" y="{top}" width="{plotw}" height="{ploth}" fill="#f8fafc" stroke="#cbd5e1"/>')
        for p in points:
            x=left+(p[0]-xmin)/(xmax-xmin)*plotw;y=top+ploth-(p[1]-ymin)/(ymax-ymin)*ploth
            body.append(f'<circle class="pt" cx="{x:.2f}" cy="{y:.2f}" r="{radius}" fill="{color(p[3],vmax)}"/>')
        for t in (0,.5,1):
            xv=xmin+t*(xmax-xmin);xp=left+t*plotw
            body.append(f'<line class="axis" x1="{xp:.1f}" y1="{top+ploth}" x2="{xp:.1f}" y2="{top+ploth+5}"/><text class="tick" x="{xp:.1f}" y="{top+ploth+18}" text-anchor="middle">{xv:.0f}</text>')
            yv=ymin+t*(ymax-ymin);yp=top+ploth-t*ploth
            body.append(f'<line class="axis" x1="{left-5}" y1="{yp:.1f}" x2="{left}" y2="{yp:.1f}"/><text class="tick" x="{left-8}" y="{yp+4:.1f}" text-anchor="end">{yv:.0f}</text>')
        body.append(f'<text class="tick" x="{left+plotw/2}" y="{top+ploth+36}" text-anchor="middle">X (source local coordinates)</text>')
        body.append(f'<text class="tick" x="{left-42}" y="{top+ploth/2}" text-anchor="middle" transform="rotate(-90 {left-42} {top+ploth/2})">Y (source local coordinates)</text>')
    # Shared color key
    lx,ly,lw,lh=310,665,580,15
    for j in range(100):
        val=vmax*j/99;body.append(f'<rect x="{lx+j*lw/100:.2f}" y="{ly}" width="{lw/100+1:.2f}" height="{lh}" fill="{color(val,vmax)}"/>')
    body.append(f'<text class="tick" x="{lx}" y="{ly+31}" text-anchor="middle">0</text>')
    for t in (.25,.5,.75,1):body.append(f'<text class="tick" x="{lx+t*lw}" y="{ly+31}" text-anchor="middle">{vmax*t:,.0f}</text>')
    body.append('<text class="tick" x="600" y="658" text-anchor="middle">Cu ppm (color position scaled with log1p)</text>')
    body.append('<text class="note" x="600" y="712" text-anchor="middle">Only 355 proxy-surface-clipped nodes are displayed. CRS, vertical datum, geological domains, and validated grid support are unavailable.</text></svg>')
    FIG.write_text("\n".join(body),encoding="utf-8")
    meta={"figure":str(FIG.relative_to(ROOT)),"observed_samples":len(sample),"surface_clipped_supported_nodes":len(estimates),"shared_color_max_ppm":vmax,"z_slice_cutoffs":cuts,"grid_node_counts_by_slice":[len(x[0]) for x in slices],"projection":"XY plan view; Z is divided into three observed-range slices for grid estimates.","limitations":["2D projection collapses points with differing Z within each plotted layer.","Only grid points that passed previous provisional support screen and proxy surface clip are displayed.","Coordinate units/CRS are not documented; the displayed grid support is not validated."]}
    META.write_text(json.dumps(meta,indent=2),encoding="utf-8")
    print(json.dumps(meta,indent=2))


if __name__=="__main__":main()
