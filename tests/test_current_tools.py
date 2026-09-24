"""Check the public interface offered by the current release."""

import json
import threading
import unittest
from http.client import HTTPConnection
from http.server import ThreadingHTTPServer
from tempfile import TemporaryDirectory
from pathlib import Path
from unittest.mock import patch

import inci_main as app


CURRENT_TOOLS = (
    'preprocessing', 'dsrna-identification', 'dsrna-plotter',
    'srna-dsrna-identification', 'srna-mapping', 'srna-control-filtering',
    'degradome-analysis', 'target-prediction', 'fasta-deduplication',
)
DEVELOPMENT_TOOLS = (
    'dsrna-adaptation', 'dsrna-enhancement', 'blast-tool',
    'rnai-susceptibility', 'orthology-inference',
)


class CurrentToolTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.workspace = TemporaryDirectory()
        cls.state_patch = patch.object(app, 'STATE_FILE', Path(cls.workspace.name) / 'state.json')
        cls.output_patch = patch.object(app, 'OUTPUT_DIR', Path(cls.workspace.name) / 'outputs')
        cls.state_patch.start()
        cls.output_patch.start()
        cls.server = ThreadingHTTPServer(('127.0.0.1', 0), app.InciRequestHandler)
        cls.worker = threading.Thread(target=cls.server.serve_forever, daemon=True)
        cls.worker.start()

    @classmethod
    def tearDownClass(cls):
        cls.server.shutdown()
        cls.server.server_close()
        cls.worker.join(timeout=2)
        cls.output_patch.stop()
        cls.state_patch.stop()
        cls.workspace.cleanup()

    def request(self, method, path, payload=None):
        connection = HTTPConnection('127.0.0.1', self.server.server_port, timeout=10)
        try:
            connection.request(method, path, body=json.dumps(payload) if payload is not None else None,
                               headers={'Content-Type': 'application/json'})
            response = connection.getresponse()
            return response.status, response.read().decode()
        finally:
            connection.close()

    def test_menu_and_each_current_tool_are_available(self):
        for path in ['/', '/settings', *('/module/' + key for key in CURRENT_TOOLS)]:
            with self.subTest(path=path):
                status, html = self.request('GET', path)
                self.assertEqual(status, 200)
                self.assertNotIn('Implementation placeholder', html)
                self.assertNotIn('This tool is not available yet.', html)
                for key in CURRENT_TOOLS:
                    self.assertIn('href="/module/' + key + '"', html)
                for key in DEVELOPMENT_TOOLS:
                    self.assertNotIn('/module/' + key, html)

    def test_development_tools_cannot_be_opened_directly(self):
        for key in DEVELOPMENT_TOOLS:
            with self.subTest(key=key):
                status, _ = self.request('GET', '/module/' + key)
                self.assertEqual(status, 404)

    def test_development_tools_cannot_start_jobs(self):
        with patch.object(app, 'start_module_job') as start_job:
            for key in DEVELOPMENT_TOOLS:
                with self.subTest(key=key):
                    status, _ = self.request('POST', '/run-module', {'module_key': key, 'params': {}})
                    self.assertEqual(status, 400)
            status, _ = self.request('POST', '/run-dsrna-adaptation', {})
            self.assertEqual(status, 404)
            start_job.assert_not_called()


if __name__ == '__main__':
    unittest.main()
