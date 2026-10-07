import { LitElement, html, nothing } from './lit.ts'
import type { TemplateResult } from './lit.ts'
import { exampleFileUrl } from './api.ts'
import type { Example, Export, Hint, Options, Output } from './api.ts'
import { formatBytes } from './traces.ts'

// What an example run includes: the result alone, with its parameters, or with both and the
// sample annotation. The server links exactly these files into the job.
// A source link named by its repository (or host) and file: "MannLabs/alphadia › output-format.md".
export function sourceLabel (url: string): string {
  const { host, pathname } = new URL(url)
  const parts = pathname.split('/').filter(Boolean)
  const site = host === 'github.com' && parts.length >= 2 ? `${parts[0]}/${parts[1]}` : host.replace(/^www\./, '')
  const file = parts.length > (host === 'github.com' ? 2 : 0) ? parts[parts.length - 1] : ''
  return file ? `${site} › ${decodeURIComponent(file)}` : site
}

export type ExampleFiles = 'data' | 'data,params' | 'data,params,annotation'

interface Variant {
  files: ExampleFiles
  label: string
}

const VARIANTS: Variant[] = [
  { files: 'data', label: 'Result file only' },
  { files: 'data,params', label: '+ parameter file' },
  { files: 'data,params,annotation', label: '+ parameters + SDRF' }
]

// Software first, then either an example or the user's own files. Emits `submit-job` with the
// multipart form the server's POST /api/jobs reads.
export class ApbUploadForm extends LitElement {
  static properties = {
    options: { attribute: false },
    examples: { attribute: false },
    busy: { type: Boolean },
    output: { state: true },
    software: { state: true },
    hasData: { state: true },
    hasParams: { state: true },
    example: { state: true },
    exampleFiles: { state: true }
  }

  options: Options | null = null
  examples: Example[] = []
  busy = false
  output = ''
  software = ''
  hasData = false
  hasParams = false
  example: Example | null = null
  exampleFiles: ExampleFiles = 'data,params,annotation'

  createRenderRoot (): this { return this }

  private chosen (): Output | undefined {
    return this.options?.outputs.find(({ id }) => id === (this.output || this.options?.outputs[0]?.id))
  }

  private hint (): Hint | undefined {
    return this.options?.hints.software[this.software]
  }

  private softwareExamples (): Example[] {
    return this.examples.filter(({ software }) => software === this.software)
  }

  private file (name: string): File | undefined {
    return this.querySelector<HTMLInputElement>(`input[name="${name}"]`)?.files?.[0]
  }

  private chooseSoftware (software: string): void {
    this.software = software
    this.example = null
    this.hasData = false
    this.hasParams = false
  }

  private unavailable (example: Example, files: ExampleFiles): string | null {
    if (files === 'data' && this.hint()?.params_required && example.label === example.software) return `${example.software} needs its parameter file`
    if (files !== 'data' && example.params === null) return 'this example has no parameter file'
    if (files === 'data,params,annotation' && example.annotation === null) return 'this example has no SDRF'
    return null
  }

  private submit (event: Event): void {
    event.preventDefault()
    const output = this.chosen()
    if (!output || !this.software) return
    const form = new FormData()
    form.append('output', output.id)
    if (this.example) {
      form.append('example', this.example.id)
      form.append('example_files', this.exampleFiles)
    } else {
      const data = [...(this.querySelector<HTMLInputElement>('input[name="data"]')?.files ?? [])]
      if (data.length === 0) return
      for (const file of data) form.append('data', file)
      form.append('software', this.software)
      const params = this.file('params')
      if (params) form.append('params', params)
      const annotation = this.file('annotation')
      if (annotation && output.annotation) form.append('annotation', annotation)
    }
    this.dispatchEvent(new CustomEvent('submit-job', { detail: form, bubbles: true }))
  }

  private renderFiles (names: string[]): TemplateResult {
    return html`${names.map((name, index) => html`${index > 0 ? ' + ' : ''}<code>${name}</code>`)}`
  }

  private renderHint (hint: Hint | undefined): TemplateResult {
    if (!hint || !this.options) return html`<p class="note">No file hint for ${this.software}.</p>`
    const required = hint.params_required === true ? 'required' : 'optional'
    return html`
      <div class="hint">
        <table>
          <thead><tr><th>Version</th><th>Result file</th><th>Parameter file · ${required}</th></tr></thead>
          <tbody>
            ${hint.exports.map(({ versions, result, params, note }) => html`
              <tr>
                <td>${versions}</td>
                <td>${this.renderFiles(result)}${note ? html`<div class="note">${note}</div>` : nothing}</td>
                <td>${params === null ? 'none' : html`<code>${params}</code>`}</td>
              </tr>
            `)}
          </tbody>
        </table>
        <p class="note">Sample annotation, optional: ${this.options.hints.annotation}</p>
        ${hint.sources?.length
          ? html`<p class="note">Names from: ${hint.sources.map((url, index) => html`${index > 0 ? ', ' : ''}<a href=${url} target="_blank" rel="noopener">${sourceLabel(url)}</a>`)}</p>`
          : nothing}
      </div>
    `
  }

  private renderExamples (examples: Example[]): TemplateResult | typeof nothing {
    if (examples.length === 0) return html`<p class="note try">No example for ${this.software} yet; upload your own files.</p>`
    return html`
      <div class="try">
        <p class="note">Try an example; each runs as-is on this server, no upload:</p>
        <table>
          ${examples.map(example => html`
            <tr>
              <td>${example.label}${example.version ? ` ${example.version}` : ''}<br><span class="note">${example.module} · ${example.folder ? `${example.data[0]}/` : example.data.join(' + ')} · ${formatBytes(example.bytes)}</span></td>
              <td><div class="variants">
                ${VARIANTS.map(({ files, label }) => {
                  const reason = this.unavailable(example, files)
                  const active = this.example === example && this.exampleFiles === files
                  return html`
                    <button type="button" class="secondary ${active ? 'active' : ''}" ?disabled=${reason !== null}
                      title=${reason ?? ''} @click=${() => { this.example = example; this.exampleFiles = files }}>
                      ${label}
                    </button>
                  `
                })}
              </div></td>
            </tr>
          `)}
        </table>
      </div>
    `
  }

  private renderExample (example: Example): TemplateResult {
    const included = new Set(this.exampleFiles.split(','))
    const storedNote = (names: string[], stored: string[]): TemplateResult | typeof nothing =>
      stored.join() !== names.join()
        ? html` <span class="note">(ProteoBench stores ${stored.length > 1 ? 'them' : 'it'} as ${stored.join(' + ')})</span>`
        : nothing
    const links = (field: string, names: string[], stored: string[]): TemplateResult | string => {
      if (names.length === 0) return '–'
      if (!included.has(field)) return html`<span class="note">not sent: ${names.join(' + ')}</span>`
      return html`${names.map((name, index) => html`${index > 0 ? ' + ' : ''}<a href=${exampleFileUrl(example.id, name)} download=${name}>${name}</a>`)}${storedNote(names, stored)}`
    }
    const optional = (value: string | null): string[] => (value === null ? [] : [value])
    const exported: Export | null = example.export
    return html`
      <div class="example">
        <dl>
          <dt>Example</dt><dd>${example.label}${example.version ? ` ${example.version}` : ''} · ${example.module}</dd>
          ${exported
            ? html`<dt>${example.software} writes</dt><dd>${this.renderFiles(exported.result)}${exported.params ? html` and <code>${exported.params}</code>` : nothing}</dd>`
            : nothing}
          <dt>Result ${example.data.length > 1 ? 'files' : 'file'}</dt>
          <dd>
            ${example.folder ? html`<code>${example.data[0]}/</code> (folder, sent whole)` : links('data', example.data, example.stored_data)}
          </dd>
          <dt>Parameter file</dt><dd>${links('params', optional(example.params), optional(example.stored_params))}</dd>
          <dt>Sample annotation</dt><dd>${links('annotation', optional(example.annotation), optional(example.annotation))}</dd>
        </dl>
        ${exported?.note ? html`<p class="note">${exported.note}</p>` : nothing}
        <button type="button" class="secondary" @click=${() => { this.example = null }}>Upload my own files instead</button>
      </div>
    `
  }

  private renderUploads (): TemplateResult {
    const output = this.chosen()
    const required = this.hint()?.params_required === true
    return html`
      <label>
        Result file <span class="optional">several when the tool writes several, e.g. AlphaDIA 1.12</span>
        <input type="file" name="data" required multiple
          @change=${(event: Event) => { this.hasData = Boolean((event.currentTarget as HTMLInputElement).files?.length) }}>
      </label>
      <label>
        Parameter file <span class="optional">${required ? 'required' : 'optional'}</span>
        <input type="file" name="params"
          @change=${(event: Event) => { this.hasParams = Boolean((event.currentTarget as HTMLInputElement).files?.length) }}>
      </label>
      <label>
        Sample annotation <span class="optional">optional</span>
        <input type="file" name="annotation" ?disabled=${output !== undefined && !output.annotation}>
      </label>
    `
  }

  private renderOutputs (chosen: Output | undefined): TemplateResult {
    const outputs = this.options?.outputs ?? []
    return html`
      <fieldset class="outputs">
        <legend>Output</legend>
        <div class="table-scroll">
          <table>
            <thead><tr><th>Format</th><th>File</th><th>What you get</th><th>Open it with</th><th>More</th></tr></thead>
            <tbody>
              ${outputs.map(({ id, about }) => html`
                <tr class=${id === chosen?.id ? 'chosen' : ''} @click=${() => { this.output = id }}>
                  <td>
                    <label class="choice">
                      <input type="radio" name="output" .checked=${id === chosen?.id} @change=${() => { this.output = id }}>
                      ${about.title}${about.version ? html` <span class="note">${about.version}</span>` : nothing}
                    </label>
                  </td>
                  <td><code>${about.file}</code></td>
                  <td>${about.what}</td>
                  <td><code>${about.open}</code></td>
                  <td>${about.links.map((link, index) => html`${index > 0 ? html`<br>` : nothing}<a href=${link.url} target="_blank" rel="noopener" @click=${(event: Event) => { event.stopPropagation() }}>${link.label}</a>`)}</td>
                </tr>
              `)}
            </tbody>
          </table>
        </div>
      </fieldset>
    `
  }

  private ready (): boolean {
    if (this.busy || !this.software) return false
    if (this.example) return true
    return this.hasData && (this.hasParams || this.hint()?.params_required !== true)
  }

  render (): TemplateResult {
    if (!this.options) return html`<p class="note">Loading options…</p>`
    const output = this.chosen()
    return html`
      <form class="card upload" @submit=${(event: Event) => { this.submit(event) }}>
        <h2>1 · Software and files</h2>
        <label>
          Software
          <select @change=${(event: Event) => { this.chooseSoftware((event.currentTarget as HTMLSelectElement).value) }}>
            <option value="" ?selected=${this.software === ''} disabled>Choose the software that wrote the result…</option>
            ${this.options.software.map(name => html`<option value=${name} ?selected=${name === this.software}>${name}</option>`)}
          </select>
        </label>
        ${this.software
          ? html`
            ${this.renderHint(this.hint())}
            ${this.renderExamples(this.softwareExamples())}
            ${this.example ? this.renderExample(this.example) : this.renderUploads()}
            ${output !== undefined && !output.annotation
              ? html`<p class="note">${output.name} takes no annotation; CV is computed across all samples.</p>`
              : nothing}
          `
          : nothing}
        ${this.renderOutputs(output)}
        <div class="actions">
          <button type="submit" ?disabled=${!this.ready()}>${this.busy ? 'Working…' : 'Convert'}</button>
          <span class="note">Uploads up to ${formatBytes(this.options.max_upload_bytes)} in total.</span>
        </div>
        <p class="note source">File names: ${this.options.hints.source}.</p>
      </form>
    `
  }
}
