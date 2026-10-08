// Composition root: define the elements, load the options, submit, poll, show results.
// The job id lives in the URL hash, so a reload resumes the same job.

import './app.css'
import { loadExamples, loadOptions, loadQc, loadStatus, submitJob } from './api.ts'
import type { Status } from './api.ts'
import { ApbAbout } from './about.ts'
import { ApbDownloads } from './downloads.ts'
import { ApbFilePreview } from './file-preview.ts'
import { ApbJobStatus } from './job-status.ts'
import { ApbQcPanel } from './qc-panel.ts'
import { ApbUploadForm } from './upload-form.ts'

customElements.define('apb-upload-form', ApbUploadForm)
customElements.define('apb-job-status', ApbJobStatus)
customElements.define('apb-downloads', ApbDownloads)
customElements.define('apb-qc-panel', ApbQcPanel)
customElements.define('apb-file-preview', ApbFilePreview)
customElements.define('apb-about', ApbAbout)

const POLL_MS = 2000

function element<T extends HTMLElement> (selector: string): T {
  const found = document.querySelector<T>(selector)
  if (!found) throw new Error(`${selector} is missing from the page`)
  return found
}

const form = element<ApbUploadForm>('apb-upload-form')
const progress = element<ApbJobStatus>('apb-job-status')
const downloads = element<ApbDownloads>('apb-downloads')
const qc = element<ApbQcPanel>('apb-qc-panel')
const about = element<ApbAbout>('apb-about')

let polling = 0

function show (status: Status | null): void {
  progress.status = status
  downloads.status = status
  form.busy = status !== null && (status.state === 'queued' || status.state === 'running')
}

async function follow (id: string): Promise<void> {
  window.clearTimeout(polling)
  progress.error = ''
  qc.qc = null
  try {
    const status = await loadStatus(id)
    show(status)
    if (status.state === 'done') qc.qc = await loadQc(id)
    else if (status.state !== 'failed') polling = window.setTimeout(() => { void follow(id) }, POLL_MS)
  } catch (error) {
    show(null)
    progress.error = String(error)
  }
}

form.addEventListener('submit-job', event => {
  const body = (event as CustomEvent<FormData>).detail
  form.busy = true
  progress.error = ''
  submitJob(body)
    .then(id => { window.location.hash = `job=${id}` })
    .catch((error: unknown) => {
      form.busy = false
      progress.error = String(error)
    })
})

function fromHash (): void {
  const id = new URLSearchParams(window.location.hash.slice(1)).get('job')
  if (id) void follow(id)
}

window.addEventListener('hashchange', fromHash)
loadOptions()
  .then(options => { form.options = options; about.versions = options.versions })
  .catch((error: unknown) => { progress.error = String(error) })
loadExamples()
  .then(examples => { form.examples = examples })
  .catch((error: unknown) => { progress.error = String(error) })
fromHash()
