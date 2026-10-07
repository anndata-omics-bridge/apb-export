import { LitElement, html, nothing } from './lit.ts'
import type { PropertyValues, TemplateResult } from './lit.ts'
import { Plotly } from './plotly.ts'
import type { Qc, QcMatrix } from './api.ts'
import { cvTraces, densityTraces, formatPercent, matrixLabel } from './traces.ts'

const ALL_SAMPLES = 'all samples'

const LAYOUT = {
  margin: { t: 40, r: 10, b: 20, l: 60 },
  paper_bgcolor: 'rgba(0,0,0,0)',
  plot_bgcolor: 'rgba(0,0,0,0)',
  font: { size: 12 },
  // Below the axis title: Plotly grows the bottom margin to fit a legend anchored at its top.
  legend: { orientation: 'h' as const, yanchor: 'top' as const, y: -0.3 }
}

// One tab per result matrix; Plotly draws only into the visible tab so it can measure it.
export class ApbQcPanel extends LitElement {
  static properties = {
    qc: { attribute: false },
    selected: { state: true }
  }

  qc: Qc | null = null
  selected = 0

  createRenderRoot (): this { return this }

  private matrix (): QcMatrix | undefined {
    return this.qc?.matrices[this.selected]
  }

  protected updated (changed: PropertyValues): void {
    if (changed.has('qc')) this.selected = Math.min(this.selected, Math.max(0, (this.qc?.matrices.length ?? 1) - 1))
    const matrix = this.matrix()
    const density = this.querySelector<HTMLElement>('[data-plot="density"]')
    const cv = this.querySelector<HTMLElement>('[data-plot="cv"]')
    if (!matrix || !density || !cv) return
    void Plotly.react(density, densityTraces(matrix), {
      ...LAYOUT,
      title: { text: 'Intensity density per sample' },
      xaxis: { title: { text: 'log2 intensity' } },
      yaxis: { title: { text: 'density' } }
    }, { responsive: true, displaylogo: false })
    void Plotly.react(cv, cvTraces(matrix), {
      ...LAYOUT,
      barmode: 'overlay',
      // Always shown: the legend carries each group's median, also for a single group.
      showlegend: true,
      title: { text: matrix.cv.grouping === ALL_SAMPLES ? 'CV per feature across all samples' : `CV per feature within each ${matrix.cv.grouping}` },
      xaxis: { title: { text: 'CV (%)' } },
      yaxis: { title: { text: 'features' } }
    }, { responsive: true, displaylogo: false })
  }

  render (): TemplateResult | typeof nothing {
    const matrices = this.qc?.matrices ?? []
    const matrix = this.matrix()
    if (!matrix) return nothing
    const empty = matrix.cv.series.every(({ features }) => features === 0)
    return html`
      <section class="card qc">
        <h2>4 · Quality control</h2>
        ${matrices.length > 1
          ? html`<nav class="tabs" role="tablist">
              ${matrices.map((item, index) => html`
                <button role="tab" aria-selected=${String(index === this.selected)} @click=${() => { this.selected = index }}>
                  ${matrixLabel(item)}
                </button>
              `)}
            </nav>`
          : html`<p class="note">${matrixLabel(matrix)}</p>`}
        <p class="summary">${matrix.features.toLocaleString()} features × ${matrix.samples.length} samples</p>
        <div class="plots">
          <div data-plot="density" class="plot"></div>
          <div data-plot="cv" class="plot"></div>
        </div>
        ${empty ? html`<p class="note">No feature is quantified in two samples of one group, so there is no CV.</p>` : nothing}
        <table class="counts">
          <thead><tr><th>Sample</th><th>Group</th><th>Detected</th><th>Missing</th></tr></thead>
          <tbody>
            ${matrix.samples.map(sample => html`
              <tr>
                <td>${sample.name}</td>
                <td>${sample.group ?? '–'}</td>
                <td>${sample.detected.toLocaleString()}</td>
                <td>${formatPercent(sample.missing_fraction)}</td>
              </tr>
            `)}
          </tbody>
        </table>
      </section>
    `
  }
}
