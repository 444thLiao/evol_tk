####################################################################################################
#
# EntrezDownloader
# Author: Leon Kuchenbecker
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
####################################################################################################

import requests
import io
import threading
import time

from concurrent.futures import ThreadPoolExecutor
from concurrent import futures

class RequestLimiter:

    def __init__(self, min_wait = 0.4):
        """The RequestLimiter class provides functionality to limit the rate at which new requests are made."""
        self.lock = threading.Lock()
        self.last_request = None
        self.min_wait = min_wait

    def wait(self):
        """The wait() function blocks until a minimum wait time from the previous invocation has passed. Thread safe."""
        with self.lock:
            # This is the first request
            if not self.last_request:
                self.last_request = time.time()
                return

            # This is not the first request
            diff = time.time() - self.last_request
            if diff < self.min_wait:
                tsleep = self.min_wait - diff
                time.sleep(tsleep)
            self.last_request = time.time()

class ResultCollector:

    def __init__(self, pbar = None):
        """The ResultCollector class provides functionality for threads to deliver their results."""
        self.pbar = pbar
        self.results = []
        self.failed = []
        self.lock = threading.Lock()

    def add_results(self, results):
        """Adds results to the collector. If a progress bar was provided, it updates the progress bar."""
        with self.lock:
            self.results += results
            if self.pbar:
                self.pbar.update(len(results))

    def add_failed(self, ids):
        """Adds failed IDs to the collector. If a progress bar was provided, it updates the progress bar."""
        with self.lock:
            self.failed += ids
            if self.pbar:
                self.pbar.update(len(ids))

class EntrezDownloader:

    def __init__(self, num_threads = 30, batch_size = 10, email = None, api_key = None, pbar = False):
        """The EntrezDownloader class enables parallel downloads via the NCBI Entrez interface"""
        self.baseurl         = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils'
        self.num_threads     = num_threads
        self.batch_size      = batch_size
        self.email           = email
        self.api_key         = api_key
        self.request_limiter = RequestLimiter(min_wait=0.35 if not api_key else 0.15)
        self.print_lock      = threading.Lock()
        self.pbar            = pbar

    def _efetch_batch(self, db, ids, result_collector, result_func, **kwargs):
        post_data = {
            'tool'    : 'EntrezDownloader',
            'email'   : self.email,
            'api_key' : self.api_key,
            'id'      : ','.join(ids),
            'db'      : db,
            'retmode' : 'xml',
        }

        post_data.update(kwargs)

        if self.email:
            post_data.update({'email':self.email})

        if self.api_key:
            post_data.update({'api_key':self.api_key})

        error = None
        for i in range(3): # Retry three times
            try:
                self.request_limiter.wait()
                response = requests.post(f'{self.baseurl}/efetch.cgi', post_data)
                if response.status_code == 200:
                    results  = result_func(response.text)
                    result_collector.add_results(results)
                    error = None
                    break
                else:
                    error = f'[STATUS {response.status_code}] An error occurred, you may see a response text here: {response.text}'
            except Exception as e:
                error = f'[UNKNOWN ERROR] {e}'

        if error:
            result_collector.add_failed(ids)
            print(error)

    def efetch(self, db, ids, result_func = lambda x : [x], **kwargs):
        """Interface to the efetch database.

        result_func: A function to be applied to the response. Must return an iterable.
        """

        if self.pbar:
            from tqdm import tqdm
            results = ResultCollector(pbar = tqdm(total=len(ids), unit = 'records'))
        else:
            results = ResultCollector()

        executor = ThreadPoolExecutor(max_workers=self.num_threads)

        fs = []
        for start in range(0, len(ids), self.batch_size):
            num = len(ids)-start
            num = self.batch_size if num > self.batch_size else num
            f = executor.submit(self._efetch_batch,
                    db = db,
                    ids = ids[start:start+num],
                    result_collector = results,
                    result_func = result_func,
                    **kwargs)
            fs.append(f)

        futures.wait(fs)

        return results.results, results.failed
