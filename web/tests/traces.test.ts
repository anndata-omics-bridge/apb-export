import assert from 'node:assert/strict'
import { test } from 'node:test'
import type { QcMatrix } from '../src/api.ts'
import { cvTraces, densityTraces, formatBytes, formatPercent, matrixLabel } from '../src/traces.ts'
import { sourceLabel } from '../src/upload-form.ts'

const matrix: QcMatrix = {
  file: 'run_prolfqua.h5ad',
  level: 'X',
  layer: 'X',
  features: 3,
  samples: [
    { name: 'a1', group: 'A', detected: 3, missing_fraction: 0 },
    { name: 'b1', group: 'B', detected: 2, missing_fraction: 1 / 3 }
  ],
  density: { x: [1, 2], series: [{ sample: 'a1', y: [0.5, 0.5] }, { sample: 'b1', y: [1, 0] }] },
  cv: {
    grouping: 'condition',
    x: [2, 6],
    series: [
      { group: 'A', samples: 2, features: 3, median: 4.25, counts: [2, 1] },
      { group: 'B', samples: 1, features: 0, median: null, counts: [0, 0] }
    ]
  }
}

test('one density line per sample on the shared bins', () => {
  const traces = densityTraces(matrix)
  assert.deepEqual(traces.map(({ name }) => name), ['a1', 'b1'])
  assert.deepEqual(traces[0]?.x, [1, 2])
})

test('CV bars only for groups with features, median in the legend', () => {
  const traces = cvTraces(matrix)
  assert.deepEqual(traces.map(({ name }) => name), ['A (median 4.3%)'])
  assert.deepEqual(traces[0]?.y, [2, 1])
})

test('labels and units', () => {
  assert.equal(matrixLabel(matrix), 'run_prolfqua.h5ad · X · X')
  assert.equal(formatBytes(512), '512 B')
  assert.equal(formatBytes(21_704_392), '21.7 MB')
  assert.equal(formatPercent(null), '–')
  assert.equal(formatPercent(0.0585), '5.9%')
})

test('source links name their repository and file', () => {
  assert.equal(sourceLabel('https://github.com/MannLabs/alphadia/blob/main/docs/methods/output-format.md'), 'MannLabs/alphadia › output-format.md')
  assert.equal(sourceLabel('https://github.com/vdemichev/DiaNN'), 'vdemichev/DiaNN')
  assert.equal(sourceLabel('https://sage-docs.vercel.app/docs/results'), 'sage-docs.vercel.app › results')
})
