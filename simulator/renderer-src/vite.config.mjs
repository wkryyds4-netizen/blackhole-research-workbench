import {defineConfig} from 'vite';
import glsl from 'vite-plugin-glsl';
export default defineConfig({plugins:[glsl()],build:{target:'es2022',outDir:'../public/lensing',emptyOutDir:true,lib:{entry:'adapter.js',formats:['es'],fileName:()=> 'renderer.js'},sourcemap:true}});
