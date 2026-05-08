import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

ADSL_PATH = 'ADSL.csv'
ADRS_PATH = 'ADRS.csv'
ADTR_PATH = 'ADTR.csv'

RESP_COLORS = {'CR':'#1f4788','PR':'#5a9bd4','SD':'#f0c929','PD':'#c41e3a'}


def build_df(adsl, adrs, adtr):
    rows = []
    for _, s in adsl.iterrows():
        uid = s['USUBJID']

        bl           = adtr[(adtr['USUBJID'] == uid) & (adtr['AVISIT'] == 'BASELINE')]
        baseline_sum = pd.to_numeric(bl['AVAL'].iloc[0], errors='coerce') if len(bl) > 0 else np.nan
        has_baseline = not pd.isna(baseline_sum) and baseline_sum > 0

        post   = adtr[(adtr['USUBJID'] == uid) & (adtr['AVISIT'] != 'BASELINE')]
        n_post = len(post)
        best_pct = np.nan
        if has_baseline and n_post > 0:
            sums = pd.to_numeric(post['AVAL'], errors='coerce').dropna()
            if len(sums) > 0:
                best_pct = ((sums.min() - baseline_sum) / baseline_sum) * 100

        pt_rs     = adrs[(adrs['USUBJID'] == uid) & (adrs['PARAMCD'] == 'OVRLRESP')]
        best_resp = 'NE'
        if len(pt_rs) > 0:
            rank  = {'CR':1,'PR':2,'SD':3,'PD':4}
            resps = [r for r in pt_rs['AVALC'].tolist() if r in rank]
            if resps:
                best_resp = min(resps, key=lambda x: rank[x])

        category = 'Evaluable' if (has_baseline and n_post >= 1) else 'Not Evaluable'
        if pd.isna(best_pct):
            best_pct = 0.0

        rows.append({
            'usubjid':    uid,
            'cohort':     s.get('COHORT', ''),
            'liver_mets': s.get('LIVERMETS', 'N') == 'Y',
            'best_resp':  best_resp,
            'best_pct':   best_pct,
            'category':   category,
        })

    return pd.DataFrame(rows)


def draw_waterfall(data, title, filename):
    data = data.sort_values('best_pct', ascending=False).reset_index(drop=True)
    n    = len(data)
    fig, ax = plt.subplots(figsize=(max(12, n * 0.5), 8))
    fig.subplots_adjust(bottom=0.22)

    y_max = max(100, data['best_pct'].max() + 15)
    y_min = min(-100, data['best_pct'].min() - 15)
    ax.set_xlim(0.5, n + 0.5)
    ax.set_ylim(y_min, y_max)
    ax.set_ylabel('Best % change from baseline in target lesion', fontsize=11, fontweight='bold')
    ax.set_title(title, fontsize=12, fontweight='bold', pad=14)

    ax.axhline(y=0,   color='black', linewidth=2,   zorder=3)
    ax.axhline(y=20,  color='black', linewidth=1.5, linestyle='--', zorder=3)
    ax.axhline(y=-30, color='black', linewidth=1.5, linestyle='--', zorder=3)
    ax.text(n - 0.5,  22,  '20% (PD threshold)', va='bottom', ha='right', fontsize=8)
    ax.text(n - 0.5, -28, '-30% (PR threshold)', va='top',   ha='right', fontsize=8)
    ax.grid(axis='y', color='gray', linestyle='--', linewidth=0.5, alpha=0.3, zorder=0)

    x_pos = np.arange(1, n + 1)
    for i, (_, row) in enumerate(data.iterrows()):
        if row['category'] == 'Not Evaluable':
            ax.bar(x_pos[i], 0.8, width=0.8, bottom=-0.4,
                   color='#555555', edgecolor='black', linewidth=1.2, zorder=2)
        else:
            ax.bar(x_pos[i], row['best_pct'], width=0.8,
                   color=RESP_COLORS.get(row['best_resp'], 'gray'),
                   edgecolor='black', linewidth=1.2, zorder=2)

    ax.set_xticks(x_pos)
    ax.set_xticklabels(data['usubjid'], rotation=90, fontsize=7)

    info_y = y_min - (y_max - y_min) * 0.08
    ax.text(0.3, info_y, 'Liver\nmets', ha='center', va='top', fontsize=8, fontweight='bold')
    for i, (_, row) in enumerate(data.iterrows()):
        ax.text(x_pos[i], info_y, 'Y' if row['liver_mets'] else 'N',
                ha='center', va='top', fontsize=7)

    evaluable = data[data['category'] == 'Evaluable']
    rc        = evaluable['best_resp'].value_counts()
    legend_elements = [
        Patch(facecolor=RESP_COLORS[r], edgecolor='black', linewidth=1.2,
              label=f'{r}  (n={rc.get(r, 0)})')
        for r in ['CR','PR','SD','PD']
    ] + [Patch(facecolor='#555555', edgecolor='black', linewidth=1.2, label='Not Evaluable')]
    legend = ax.legend(handles=legend_elements, loc='upper right',
                       title='Best Response', title_fontsize=10,
                       fontsize=9, frameon=True, fancybox=False, edgecolor='black')
    legend.get_frame().set_linewidth(1.5)

    plt.savefig(filename, dpi=300, bbox_inches='tight')
    plt.close()
    print(f'Saved: {filename}')


def main():
    adsl = pd.read_csv(ADSL_PATH, parse_dates=['TRTSDT','TRTEDT','CUTDTC'])
    adrs = pd.read_csv(ADRS_PATH, parse_dates=['ADT'])
    adtr = pd.read_csv(ADTR_PATH, parse_dates=['ADT'])

    df        = build_df(adsl, adrs, adtr)
    evaluable = df[df['category'] == 'Evaluable']
    rc        = evaluable['best_resp'].value_counts()

    draw_waterfall(
        evaluable,
        f'Waterfall Plot — Evaluable Patients (N={len(evaluable)})\n'
        f'CR:{rc.get("CR",0)}  PR:{rc.get("PR",0)}  SD:{rc.get("SD",0)}  PD:{rc.get("PD",0)}',
        'waterfall_evaluable.png'
    )

    draw_waterfall(df, f'Waterfall Plot — All Patients (N={len(df)})', 'waterfall_all.png')

    for cohort in sorted(df['cohort'].dropna().unique()):
        sub  = df[df['cohort'] == cohort]
        ev   = sub[sub['category'] == 'Evaluable']
        rc_c = ev['best_resp'].value_counts()
        draw_waterfall(
            sub,
            f'Waterfall Plot — Cohort {int(cohort)} (N={len(sub)})\n'
            f'CR:{rc_c.get("CR",0)}  PR:{rc_c.get("PR",0)}  SD:{rc_c.get("SD",0)}  PD:{rc_c.get("PD",0)}',
            f'waterfall_cohort{int(cohort)}.png'
        )


if __name__ == '__main__':
    main()
