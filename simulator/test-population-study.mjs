// Offline release checks for the saved study and its display contract.
import assert from 'node:assert/strict';
import {readFileSync} from 'node:fs';
import {createHash} from 'node:crypto';
import {fileURLToPath} from 'node:url';
import {validatePopulationStudy, selectedStudyResult, populationStudySession, medianEnvelope} from './public/population-study-contract.mjs';

const root = new URL('../', import.meta.url);
const read = path => readFileSync(new URL(path, root));
const sha = bytes => createHash('sha256').update(bytes).digest('hex');
const bytes = read('simulator/public/datasets/population-study.json');
const data = validatePopulationStudy(JSON.parse(bytes));
const manifest = JSON.parse(read('research/population-invariant/workbench-export.json'));
assert.equal(sha(bytes), manifest.sha256);
for (const report of manifest.source_reports) assert.equal(sha(read(report.repo_path)), report.sha256, report.repo_path);
for (const paper of manifest.papers) assert.equal(sha(read('simulator/public/' + paper.url)), paper.sha256, paper.url);
assert.equal(sha(read('research/population-invariant/export_workbench.py')), manifest.exporter_sha256);
assert.equal(data.observations.events, 69);
assert.equal(data.observations.posterior_rows, 11180578);
assert.equal(data.paired_banks.replicates, 512);
assert.equal(selectedStudyResult(data, data.defaults).assumed, null);

let selections = 0;
for (const metric of data.metrics) for (const model of metric.models) for (const policy of data.policies) {
  const settings = {metric_id: metric.id, model_id: model.id, policy_id: policy.id, assumed_B: 0.5};
  const result = selectedStudyResult(data, settings);
  assert.equal(result.contrast.metric_id, metric.id);
  assert.equal(result.assumed.estimated_from_data, false);
  assert.ok(result.assumed.median_bounds[0] <= result.posterior.quantiles[1]);
  assert.ok(result.assumed.median_bounds[1] >= result.posterior.quantiles[1]);
  const session = JSON.parse(JSON.stringify(populationStudySession(data, settings, sha(bytes))));
  validatePopulationStudy(session.study);
  assert.deepEqual(selectedStudyResult(session.study, session.settings).posterior.quantiles, result.posterior.quantiles);
  selections++;
}
assert.equal(selections, 24);
// Analytic uniform-density reference, independent of the stored model curves.
const envelope = medianEnvelope([0, 0.5, 1], [1, 1, 1], Math.log(3));
assert.ok(Math.abs(envelope[0] - 0.25) < 1e-14);
assert.ok(Math.abs(envelope[1] - 0.75) < 1e-14);

const reject = mutate => { const copy = structuredClone(data); mutate(copy); assert.throws(() => validatePopulationStudy(copy)); };
reject(copy => { copy.support.estimated_from_data = true; });
reject(copy => { copy.paired_banks.contrasts[0].metric_id = 'unknown'; });
reject(copy => { copy.metrics[0].models[0].results[data.policies[0].id].quantiles[1] += 0.01; });
reject(copy => { copy.support.tail_plan.executed = true; });
console.log(`PASS: population release hashes, ${selections} selections/session round-trips, analytic envelope and invalid scientific labels (${fileURLToPath(root)}).`);
