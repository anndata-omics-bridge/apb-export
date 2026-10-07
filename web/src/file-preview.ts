import { LitElement, html, nothing } from './lit.ts'
import type { TemplateResult } from './lit.ts'
import { exampleFileUrl, loadHead } from './api.ts'
import type { Head } from './api.ts'
import { formatBytes } from './traces.ts'

// A modal showing the first lines of an example file, so a visitor sees what the file looks
// like without downloading it; the download sits below.
export class ApbFilePreview extends LitElement {
  static properties = {
    head: { state: true },
    error: { state: true }
  }

  head: Head | null = null
  error = ''
  private url = ''

  createRenderRoot (): this { return this }

  async show (exampleId: string, name: string): Promise<void> {
    this.url = exampleFileUrl(exampleId, name)
    this.head = null
    this.error = ''
    this.querySelector('dialog')?.showModal()
    try {
      this.head = await loadHead(exampleId, name)
    } catch (error) {
      this.error = String(error)
    }
  }

  private close (): void {
    this.querySelector('dialog')?.close()
  }

  private renderHead (head: Head): TemplateResult {
    const shown = head.format === 'parquet'
      ? `first ${head.lines.length - 1} rows of a Parquet file, shown as tab-separated text`
      : head.complete ? 'the whole file' : `first ${head.lines.length} lines`
    return html`
      <p class="note">${formatBytes(head.bytes)} · ${shown}</p>
      <pre class="head">${head.lines.join('\n')}</pre>
    `
  }

  render (): TemplateResult {
    const name = this.head?.name ?? decodeURIComponent(this.url.split('/').pop() ?? '')
    return html`
      <dialog class="preview" @click=${(event: Event) => { if (event.target === event.currentTarget) this.close() }}>
        <div class="preview-body">
          <h3><code>${name}</code></h3>
          ${this.head ? this.renderHead(this.head) : nothing}
          ${!this.head && !this.error ? html`<p class="note">Loading…</p>` : nothing}
          ${this.error ? html`<pre class="error">${this.error}</pre>` : nothing}
          <div class="preview-actions">
            <a class="button" href=${this.url} download=${name}>Download ${name}${this.head ? ` (${formatBytes(this.head.bytes)})` : ''}</a>
            <button type="button" class="secondary" @click=${() => { this.close() }}>Close</button>
          </div>
        </div>
      </dialog>
    `
  }
}
