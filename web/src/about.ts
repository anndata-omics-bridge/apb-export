import { LitElement, html, nothing } from './lit.ts'
import type { TemplateResult } from './lit.ts'

const ORGANIZATION = 'https://github.com/anndata-omics-bridge'
const APB2_DOCS = 'https://anndata-omics-bridge.github.io/apb2/'
const APB2_SOURCE = 'https://github.com/anndata-omics-bridge/apb2'
const PROTEOBENCH = 'https://github.com/Proteobench/ProteoBench'

// The header's About link and the dialog it opens: what the page is, where its parts live,
// where the examples come from, and the installed versions.
export class ApbAbout extends LitElement {
  static properties = {
    versions: { attribute: false }
  }

  versions: Record<string, string> = {}

  createRenderRoot (): this { return this }

  private open (): void {
    this.querySelector('dialog')?.showModal()
  }

  private close (): void {
    this.querySelector('dialog')?.close()
  }

  private link (url: string, label: string): TemplateResult {
    return html`<a href=${url} target="_blank" rel="noopener">${label}</a>`
  }

  render (): TemplateResult {
    const versions = Object.entries(this.versions)
    return html`
      <button type="button" class="secondary" @click=${() => { this.open() }}>About</button>
      <dialog class="preview" @click=${(event: Event) => { if (event.target === event.currentTarget) this.close() }}>
        <div class="about-body">
          <h3>About APB Export</h3>
          <p>
            Converts one vendor result into an APB2 result or into the input a downstream tool reads,
            and shows quality control before the download. Conversion is apb2, the AnnData Proteomics
            Bridge: a rule document per tool and version declares how its tables become AnnData
            levels.
          </p>
          <ul>
            <li>${this.link(ORGANIZATION, 'anndata-omics-bridge')} on GitHub</li>
            <li>apb2: ${this.link(APB2_DOCS, 'documentation')} and ${this.link(APB2_SOURCE, 'source')}</li>
          </ul>
          <h4>ProteoBench</h4>
          <p>
            Most examples are public ${this.link(PROTEOBENCH, 'ProteoBench')} submissions, shown under the
            names the tools write rather than ProteoBench's storage names, and annotated with the SDRFs of
            ProteoBench's modules; the pb_custom and MetaMorpheus examples are apb2's test samples.
            apb2's vendor parsing started from ProteoBench's.
          </p>
          <h4>Your files</h4>
          <p>Uploads are converted on this server and deleted with their results after a day.</p>
          ${versions.length > 0
            ? html`<p class="note">${versions.map(([name, version], index) => html`${index > 0 ? ' · ' : ''}${name} ${version}`)}</p>`
            : nothing}
          <div class="preview-actions">
            <button type="button" class="secondary" @click=${() => { this.close() }}>Close</button>
          </div>
        </div>
      </dialog>
    `
  }
}
