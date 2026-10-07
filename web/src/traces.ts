// Pure projections from qc.json to Plotly traces and display strings; no DOM.

import type { QcMatrix } from './api.ts'

export interface LineTrace {
  type: 'scatter'
  mode: 'lines'
  name: string
  x: number[]
  y: number[]
  line: { width: number }
}

export interface BarTrace {
  type: 'bar'
  name: string
  x: number[]
  y: number[]
  opacity: number
}

export function densityTraces (matrix: QcMatrix): LineTrace[] {
  return matrix.density.series.map(({ sample, y }) => ({
    type: 'scatter',
    mode: 'lines',
    name: sample,
    x: matrix.density.x,
    y,
    line: { width: 1.5 }
  }))
}

export function cvTraces (matrix: QcMatrix): BarTrace[] {
  return matrix.cv.series
    .filter(({ features }) => features > 0)
    .map(({ group, median, counts }) => ({
      type: 'bar',
      name: median === null ? group : `${group} (median ${median.toFixed(1)}%)`,
      x: matrix.cv.x,
      y: counts,
      opacity: matrix.cv.series.length > 1 ? 0.55 : 0.9
    }))
}

export function matrixLabel (matrix: QcMatrix): string {
  return `${matrix.file} · ${matrix.level} · ${matrix.layer}`
}

export function formatBytes (bytes: number): string {
  const units = ['B', 'kB', 'MB', 'GB']
  let value = bytes
  let unit = 0
  while (value >= 1000 && unit < units.length - 1) {
    value /= 1000
    unit += 1
  }
  return `${unit === 0 ? value : value.toFixed(1)} ${units[unit]}`
}

export function formatPercent (fraction: number | null): string {
  return fraction === null ? '–' : `${(100 * fraction).toFixed(1)}%`
}
