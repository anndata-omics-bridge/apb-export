// The server's files and the one upload route. Every response is a file a job or the
// server's startup wrote; nothing here asks the server to compute.

export interface OutputLink {
  label: string
  url: string
}

export interface About {
  group: 'tool' | 'apb2'
  title: string
  file: string
  what: string
  open: string
  version: string | null
  links: OutputLink[]
}

export interface Output {
  id: string
  label: string
  command: 'apb-export' | 'apb2'
  name: string
  extension: string
  annotation: boolean
  about: About
}

export interface Export {
  rules: string[]
  versions: string
  result: string[]
  params: string | null
  kinds: string[]
  labels: string[]
  note: string | null
}

export interface Hint {
  exports: Export[]
  params_required?: boolean
  sources?: string[]
}

export interface Options {
  software: string[]
  hints: { source: string, annotation: string, software: Record<string, Hint> }
  outputs: Output[]
  max_upload_bytes: number
}

export interface Example {
  id: string
  label: string
  software: string
  version: string | null
  module: string
  data: string[]
  stored_data: string[]
  folder: boolean
  // The files inside a folder result, each previewable; empty for a file result.
  folder_files: string[]
  params: string | null
  stored_params: string | null
  annotation: string | null
  export: Export | null
  bytes: number
}

export type State = 'queued' | 'running' | 'done' | 'failed'

export interface Step {
  name: string
  state: State
  seconds: number | null
}

export interface DownloadFile {
  name: string
  bytes: number
}

export interface Status {
  id: string
  state: State
  steps: Step[]
  error: string | null
  files: DownloadFile[]
}

export interface SampleCount {
  name: string
  group: string | null
  detected: number
  missing_fraction: number | null
}

export interface CvSeries {
  group: string
  samples: number
  features: number
  median: number | null
  counts: number[]
}

export interface QcMatrix {
  file: string
  level: string
  layer: string
  features: number
  samples: SampleCount[]
  density: { x: number[], series: Array<{ sample: string, y: number[] }> }
  cv: { grouping: string, x: number[], series: CvSeries[] }
}

export interface Qc {
  matrices: QcMatrix[]
}

// The first lines of an example file, written by the server at startup; a Parquet file's
// first rows come as tab-separated text.
export interface Head {
  name: string
  bytes: number
  format: 'text' | 'parquet'
  lines: string[]
  complete: boolean
}

async function json<T> (url: string): Promise<T> {
  const response = await fetch(url, { cache: 'no-store' })
  if (!response.ok) throw new Error(`${url}: ${response.status} ${response.statusText}`)
  return await response.json() as T
}

export const loadOptions = async (): Promise<Options> => await json<Options>('api/options.json')
export const loadExamples = async (): Promise<Example[]> =>
  (await json<{ examples: Example[] }>('api/examples.json')).examples
export const exampleFileUrl = (id: string, name: string): string =>
  `api/examples/${id}/${encodeURIComponent(name)}`
export const loadHead = async (id: string, name: string): Promise<Head> =>
  await json<Head>(`api/examples/${id}/head/${encodeURIComponent(name)}`)
export const loadStatus = async (id: string): Promise<Status> => await json<Status>(`api/jobs/${id}/status.json`)
export const loadQc = async (id: string): Promise<Qc> => await json<Qc>(`api/jobs/${id}/qc.json`)
export const logUrl = (id: string): string => `api/jobs/${id}/job.log`
export const fileUrl = (id: string, name: string): string =>
  `api/jobs/${id}/files/${encodeURIComponent(name)}`

export async function submitJob (form: FormData): Promise<string> {
  const response = await fetch('api/jobs', { method: 'POST', body: form })
  const body = await response.json() as { id?: string, detail?: unknown }
  if (!response.ok || body.id === undefined) {
    const detail = typeof body.detail === 'string' ? body.detail : JSON.stringify(body.detail)
    throw new Error(`upload refused: ${detail}`)
  }
  return body.id
}
