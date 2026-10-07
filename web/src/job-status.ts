import { LitElement, html, nothing } from './lit.ts'
import type { TemplateResult } from './lit.ts'
import { logUrl } from './api.ts'
import type { Status } from './api.ts'

export class ApbJobStatus extends LitElement {
  static properties = {
    status: { attribute: false },
    error: { type: String }
  }

  status: Status | null = null
  error = ''

  createRenderRoot (): this { return this }

  render (): TemplateResult | typeof nothing {
    if (!this.status && !this.error) return nothing
    const status = this.status
    return html`
      <section class="card">
        <h2>Progress</h2>
        ${status
          ? html`
            <p class="state state-${status.state}">${status.state}</p>
            <ol class="steps">
              ${status.steps.map(step => html`
                <li class="state-${step.state}">
                  ${step.name}<span>${step.state}${step.seconds === null ? '' : ` · ${step.seconds} s`}</span>
                </li>
              `)}
            </ol>
            ${status.error ? html`<pre class="error">${status.error}</pre>` : nothing}
            <a href=${logUrl(status.id)} target="_blank" rel="noopener">Job log</a>
          `
          : nothing}
        ${this.error ? html`<pre class="error">${this.error}</pre>` : nothing}
      </section>
    `
  }
}
