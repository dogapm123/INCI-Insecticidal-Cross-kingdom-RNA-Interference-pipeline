const assert = require('node:assert/strict');
const vm = require('node:vm');
const {execFileSync} = require('node:child_process');

const source = execFileSync('python3', ['-c',
    'import inci_main; print(inci_main.render_scripts())'], {encoding: 'utf8'})
    .replace(/^[\s\S]*?<script>/, '').replace(/<\/script>[\s\S]*$/, '');

async function check(status, responseStatus = 200) {
    const values = new Map();
    const timers = [];
    const messages = [];
    const context = vm.createContext({
        document: {body: {dataset: {project: '/project'}}, querySelectorAll: () => [], addEventListener: () => {}},
        window: {location: {pathname: '/module/srna-mapping'}, setTimeout: (fn) => timers.push(fn)},
        sessionStorage: {setItem: (k, v) => values.set(k, v), getItem: k => values.get(k), removeItem: k => values.delete(k)},
        AbortSignal, Date, URL,
        fetch: async () => ({ok: responseStatus === 200, status: responseStatus,
            json: async () => ({ok: responseStatus === 200, job: {status, label: 'Mapping', logs: [], result: {message: 'Failure details'}}})})
    });
    vm.runInContext(source, context);
    context.setStatus = (message, level) => messages.push({message, level});
    context.writeTerminal = () => {};
    await context.pollPipelineJob('example');
    const key = 'inci-active-job:/project:/module/srna-mapping';
    if (status === 'running' && responseStatus === 200) {
        assert(values.has(key));
        assert.equal(timers.length, 1);
    } else {
        assert(!values.has(key));
        assert.equal(messages.at(-1).level, status === 'finished' && responseStatus === 200 ? 'success' : 'error');
        assert.equal(timers.length, status === 'finished' && responseStatus === 200 ? 1 : 0);
    }
}

(async () => {
    for (const status of ['running', 'finished', 'failed', 'stopped']) await check(status);
    await check('running', 404);
    console.log('Polling checks passed: running, completed, failed, stopped, server restart.');
})().catch(error => { console.error(error); process.exitCode = 1; });
