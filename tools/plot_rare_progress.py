"""Plot a frozen `rare_events.monitor` JSON using real TeX and observed data only.
Dependencies: matplotlib, PyMuPDF; TeX must be on PATH. No network or HPC actions.
"""
from pathlib import Path
import argparse
import hashlib
import json
import math
import shutil
import os


def plot(source,output):
    output=Path(output);output.mkdir(parents=True,exist_ok=True)
    os.environ.setdefault('MPLCONFIGDIR',str(output/'mpl_cache'))
    assert shutil.which('latex') and shutil.which('dvipng'),'A real TeX installation is required'
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import fitz
    def apply_academic_style():
        plt.rcParams.update({'text.usetex':True,'font.family':'serif','font.size':20,
            'axes.labelsize':20,'axes.titlesize':20,'xtick.labelsize':20,'ytick.labelsize':20,
            'legend.fontsize':20,'figure.titlesize':20,'axes.grid':False,
            'axes.spines.top':False,'axes.spines.right':False,'xtick.direction':'out','ytick.direction':'out'})
    apply_academic_style()
    s=json.loads(Path(source).read_text(encoding='utf-8-sig'))
    modes=[m for m in ('ams','dmc','qdmc','direct') if m in s['runs']]
    colors=dict(ams='#000000',dmc='#D62728',qdmc='#1F77B4',direct='#666666')
    styles=dict(ams='-',dmc='--',qdmc=':',direct='-.')
    def axes(title,xlabel,ylabel):
        fig,ax=plt.subplots(figsize=(8,6));fig.subplots_adjust(left=.18,right=.95,bottom=.17,top=.69)
        fig.suptitle(title,y=.97);ax.set(xlabel=xlabel,ylabel=ylabel);ax.grid(False,which='both')
        return fig,ax
    def legend(fig):
        from matplotlib.lines import Line2D
        fig.legend(handles=[Line2D([],[],color=colors[m],ls=styles[m],label=m.upper()) for m in modes],
                   loc='upper center',bbox_to_anchor=(.55,.93),ncol=2,frameon=False)
    def export(fig,name):
        fig.canvas.draw();renderer=fig.canvas.get_renderer()
        for t in fig.findobj(matplotlib.text.Text):
            if t.get_visible() and t.get_text():
                b=t.get_window_extent(renderer)
                assert b.x0>=-1 and b.y0>=-1 and b.x1<=fig.bbox.x1+1 and b.y1<=fig.bbox.y1+1,t.get_text()
        fig.savefig(output/(name+'.png'),dpi=600);fig.savefig(output/(name+'.pdf'))
        with fitz.open(output/(name+'.pdf')) as d:d[0].get_pixmap(matrix=fitz.Matrix(1.5,1.5)).save(output/(name+'_pdf_check.png'))
        plt.close(fig)
    fig,ax=axes('Rupture-reference pilot',r'Local code time $t$',r'$h_{\min}$');maxt=0.;points=0
    for m in modes:
        low={}
        for c in s['runs'][m]['calls']:
            rows=c['observations']
            if not rows:continue
            ax.plot([r[1] for r in rows],[r[2] for r in rows],color=colors[m],ls=styles[m],alpha=.15,lw=.6)
            for r in rows:low[r[1]]=min(low.get(r[1],1.),r[2]);maxt=max(maxt,r[1]);points+=1
        if low:ax.plot(sorted(low),[low[t] for t in sorted(low)],color=colors[m],ls=styles[m],lw=1.5)
    assert points>0,'Wait for real observations; no future or fabricated curves'
    ax.axhline(.01,color='black',ls=':',lw=.7);ax.set(xlim=(0,maxt*1.02),ylim=(-.02,1.03))
    legend(fig);export(fig,'hmin')
    fig,ax=axes('Compensated particle contributions',r'Selection index $k$',r'$\log_{10}(Nw_i)$');last=0
    for m in modes:
        entries=s['runs'][m]['selections'];x=[0];means=[0.];lo=[0.];hi=[0.]
        for e in entries:
            k=e.get('stage',e.get('round'))
            if m=='ams':values=[e['log_weight']/math.log(10)]
            else:values=[(e['statistics']['log_normalizer']-e['new_V'][p])/math.log(10) for p in e['parents']]
            x.append(k);means.append(sum(values)/len(values));lo.append(min(values));hi.append(max(values));last=max(last,k)
        ax.fill_between(x,lo,hi,color=colors[m],alpha=.12)
        ax.plot(x,means,color=colors[m],ls=styles[m],marker='o',fillstyle='none')
    ax.set_xlim(-.2,max(1,last*1.03))
    if not last:ax.set_ylim(-.1,.1)
    legend(fig);export(fig,'weights')
    fig,ax=axes('Measured and active solver cost','Method','Core hours');done=[];active=[]
    for m in modes:
        calls=s['runs'][m]['calls']
        done.append(sum(c['execution']['wall_seconds']*c['request']['ranks']/3600 for c in calls if c['execution']))
        active.append(sum(c['active_elapsed_core_hours'] for c in calls if not c['execution']))
    ax.bar(range(len(modes)),done,color=[colors[m] for m in modes],width=.55)
    ax.bar(range(len(modes)),active,bottom=done,color='none',edgecolor=[colors[m] for m in modes],hatch='///',width=.55)
    ax.set_xticks(range(len(modes)),[m.upper() for m in modes]);ax.set_ylim(0,max(1,max(a+b for a,b in zip(done,active))*1.3))
    export(fig,'cost')
    report=dict(job=s['job'],snapshot_utc=s['utc'],snapshot_sha256=hashlib.sha256(Path(source).read_bytes()).hexdigest(),
        observed_points=points,last_local_time=maxt,completed_solver_core_hours=dict(zip(modes,done)),
        active_elapsed_core_hours_estimate=dict(zip(modes,active)),complete_MFPT=None,
        captions=dict(hmin='bold sampled minimum envelope, not one trajectory or mean',
                      weights='range and log-weight mean, not a confidence interval; AMS provisional W is not P',
                      cost='solid completed calls; hatched elapsed active estimates including startup; no efficiency ranking'))
    (output/'plot_verification.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--snapshot',required=True);p.add_argument('--output',required=True)
    a=p.parse_args();plot(a.snapshot,a.output)
