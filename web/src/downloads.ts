import { LitElement, html, nothing } from './lit.ts'
import type { TemplateResult } from './lit.ts'
import { fileUrl } from './api.ts'
import type { Status } from './api.ts'
import { formatBytes } from './traces.ts'

export class ApbDownloads extends LitElement {
  static properties = {
    status: { attribute: false }
  }

  status: Status | null = null

  createRenderRoot (): this { return this }

  render (): TemplateResult | typeof nothing {
    const status = this.status
    if (!status || status.state !== 'done') return nothing
    return html`
      <section class="card">
        <h2>3 · Download</h2>
        <ul class="files">
          ${status.files.map(({ name, bytes }) => html`
            <li><a href=${fileUrl(status.id, name)} download=${name}>${name}</a><span>${formatBytes(bytes)}</span></li>
          `)}
        </ul>
        <p class="note">Files are deleted from the server after a day.</p>
      </section>
    `
  }
}
