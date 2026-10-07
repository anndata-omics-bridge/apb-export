import { defineConfig } from 'vite'
import { resolve } from 'node:path'

// The bundle is package data: apb-export-web serves it, so no Node runs in production.
export default defineConfig({
  base: './',
  build: {
    outDir: resolve(import.meta.dirname, '../src/apb_export/web/static'),
    emptyOutDir: true,
    assetsInlineLimit: 0
  },
  server: {
    proxy: {
      '/api': 'http://127.0.0.1:8770'
    }
  }
})
