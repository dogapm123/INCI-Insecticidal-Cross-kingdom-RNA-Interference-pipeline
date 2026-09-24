import threading
import unittest
from unittest.mock import patch

import inci_main as app


class InterfaceJobTests(unittest.TestCase):
    def setUp(self):
        self.jobs = patch.dict(app.PIPELINE_JOBS, {}, clear=True)
        self.jobs.start()
        self.addCleanup(self.jobs.stop)

    def test_abandoned_worker_is_reported_as_failed(self):
        job_id = app.create_pipeline_job('Test analysis')
        snapshot = app.pipeline_job_snapshot(job_id)
        self.assertEqual(snapshot['status'], 'failed')
        self.assertIn('stopped unexpectedly', snapshot['result']['message'])

    def test_running_worker_is_not_mistaken_for_failure(self):
        job_id = app.create_pipeline_job('Test analysis')
        release = threading.Event()
        worker = threading.Thread(target=release.wait, name='inci-test-' + job_id[:8])
        worker.start()
        try:
            self.assertEqual(app.pipeline_job_snapshot(job_id)['status'], 'running')
        finally:
            release.set()
            worker.join()

    def test_cancelled_job_cannot_report_success(self):
        job_id = app.create_pipeline_job('Test analysis')
        app.PIPELINE_JOBS[job_id]['cancel_requested'] = True
        app.finish_pipeline_job(job_id, 'finished', {}, 'Done')
        self.assertEqual(app.pipeline_job_snapshot(job_id)['status'], 'stopped')

    def test_module_failure_reaches_terminal_status(self):
        with patch.object(app, 'load_state', return_value={}), patch.object(
            app, 'run_module', side_effect=RuntimeError('Example failure')
        ):
            job_id = app.start_module_job('fasta-deduplication', {})
            for worker in threading.enumerate():
                if worker.name.endswith(job_id[:8]):
                    worker.join(timeout=2)
            snapshot = app.pipeline_job_snapshot(job_id)
            self.assertEqual(snapshot['status'], 'failed')
            self.assertEqual(snapshot['result']['message'], 'Example failure')


if __name__ == '__main__':
    unittest.main()
