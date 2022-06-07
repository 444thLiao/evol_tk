# -*- coding: utf-8 -*-
"modified by thliao at 20191019"

import requests
import io
import threading
import time

from concurrent.futures import ThreadPoolExecutor
from concurrent import futures
import pandas as pd
from bs4 import BeautifulSoup
from tqdm import tqdm


def check_retmode(retmode):
    if retmode not in ['xml','text','asn.1']:
        print("The retmode you pass might be error one. check it ")
        
class RequestLimiter:
    def __init__(self, min_wait=0.4):
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
    def __init__(self, pbar=None):
        """The ResultCollector class provides functionality for threads to deliver their results."""
        self.pbar = pbar
        self.results = []
        self.failed = []
        self.lock = threading.Lock()

    def add_results(self, results,num=None):
        """Adds results to the collector. If a progress bar was provided, it updates the progress bar."""
        with self.lock:
            self.results += results
            if self.pbar:
                if num is None:
                    self.pbar.update(len(results))
                else:
                    self.pbar.update(num)
    def add_failed(self, ids):
        """Adds failed IDs to the collector. If a progress bar was provided, it updates the progress bar."""
        with self.lock:
            self.failed += ids
            if self.pbar:
                self.pbar.update(len(ids))


class EntrezDownloader:
    def __init__(
        self, num_threads=30, batch_size=10, email=None, api_key=None, pbar=False
    ):
        """The EntrezDownloader class enables parallel downloads via the NCBI Entrez interface"""
        self.baseurl = r"https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
        self.num_threads = num_threads
        self.batch_size = batch_size
        self.email = email
        self.api_key = api_key
        self.request_limiter = RequestLimiter(min_wait=0.4 if not api_key else 0.2)
        self.print_lock = threading.Lock()
        self.pbar = pbar
        
    def disable_pbar(self):
        self.pbar = False

    def _init_pb(self, ids, batch_size=1):
        tqdm.write('it might over the progress bar since it also retrieve information of identical locus.')
        if isinstance(ids, str) and "," in ids:
            ids = [_.strip() for _ in ids.split(",") if _]
        else:
            ids = ids
        if batch_size == 1:
            batch_size = self.batch_size
        
        if self.pbar:
            results = ResultCollector(
                pbar=tqdm(total=(len(ids)) , unit="records")
            )
        else:
            results = ResultCollector()
        return ids, results
    
    def _general_batch(self, db, ids, result_collector, result_func, emode, **kwargs):

        post_data = {"email": self.email, "api_key": self.api_key, "db": db}
        if emode == "esummary":
            if isinstance(ids, str):
                ids = [ids]
            post_data["id"] = ",".join(ids)
        elif emode == "esearch":
            post_data["term"] = ids

        post_data.update(kwargs)

        if self.email:
            post_data.update({"email": self.email})

        if self.api_key:
            post_data.update({"api_key": self.api_key})

        error = None
        for i in range(3):  # Retry three times
            try:
                self.request_limiter.wait()
                response = requests.post(f"{self.baseurl}/{emode}.fcgi", post_data)
                if response.status_code == 200:
                    results = result_func(response.text)
                    if emode == "esearch" and " OR " not in ids:
                        result_collector.add_results(
                            list(zip([_.strip() for _ in ids.split(" OR ")], results))
                        )
                    else:
                        result_collector.add_results(results)
                    error = None
                    break
                else:
                    error = f"[STATUS {response.status_code}] An error occurred, you may see a response text here: {response.text}"
            except Exception as e:
                error = f"[UNKNOWN ERROR] {e}"
            except KeyboardInterrupt:
                break

        if error:
            result_collector.add_failed(ids)
            print(error)

    def _efetch_batch(
        self, db, ids, result_collector, result_func, retmode, retype, **kwargs    ):
        
        if not retmode:
            retmode = "text"
        if not retype:
            retype = "gb"
        post_data = {
            "email": self.email,
            "api_key": self.api_key,
            "id": ",".join(list(map(str, ids))),
            "db": db,
            "retmode": retmode,
            "rettype": retype,
        }
        check_retmode(retmode)
        post_data.update(kwargs)
        if self.email:
            post_data.update({"email": self.email})

        if self.api_key:
            post_data.update({"api_key": self.api_key})
        error = None
        for i in range(3):  # Retry three times
            try:
                self.request_limiter.wait()
                response = requests.post(f"{self.baseurl}/efetch.fcgi", post_data)
                if response.status_code == 200:
                    results = result_func(response.text)
                    num = len(post_data['id'].split(','))
                    result_collector.add_results(results,num=num)
                    error = None
                    break
                else:
                    error = f"[STATUS {response.status_code}] An error occurred, you may see a response text here: {response.text}"
            except Exception as e:
                error = f"[UNKNOWN ERROR] {e}"
            except KeyboardInterrupt:
                break

        if error:
            result_collector.add_failed(ids)
            print(error)

    def _esummary_batch(self, db, ids, result_collector, result_func, **kwargs):
        post_data = {
            "email": self.email,
            "api_key": self.api_key,
            "id": ",".join(ids),
            "db": db,
        }

        post_data.update(kwargs)

        if self.email:
            post_data.update({"email": self.email})

        if self.api_key:
            post_data.update({"api_key": self.api_key})

        error = None
        for i in range(3):  # Retry three times
            try:
                self.request_limiter.wait()
                response = requests.post(f"{self.baseurl}/esummary.fcgi", post_data)
                if response.status_code == 200:
                    results = result_func(response.text)
                    result_collector.add_results(results)
                    error = None
                    break
                else:
                    error = f"[STATUS {response.status_code}] An error occurred, you may see a response text here: {response.text}"
            except Exception as e:
                error = f"[UNKNOWN ERROR] {e}"
            except KeyboardInterrupt:
                break
        if error:
            result_collector.add_failed(ids)
            print(error)

    def _elink_batch(self, dbfrom, db, ids, result_collector, result_func, **kwargs):
        post_data = {
            "email": self.email,
            "api_key": self.api_key,
            "id": ids,
            "db": db,
            "dbfrom": dbfrom,
        }

        post_data.update(kwargs)

        if self.email:
            post_data.update({"email": self.email})

        if self.api_key:
            post_data.update({"api_key": self.api_key})

        error = None
        for i in range(3):  # Retry three times
            try:
                self.request_limiter.wait()
                response = requests.post(f"{self.baseurl}/elink.fcgi", post_data)
                if response.status_code == 200:
                    results = result_func(response.text)
                    result_collector.add_results(results)
                    error = None
                    break
                else:
                    error = f"[STATUS {response.status_code}] An error occurred, you may see a response text here: {response.text}"
            except Exception as e:
                error = f"[UNKNOWN ERROR] {e}"
            except KeyboardInterrupt:
                break
        if error:
            result_collector.add_failed(ids)
            print(error)


    def elink(self, dbfrom, db, ids, batch_size=1, result_func=lambda x: [x], **kwargs):
        """Interface to the elink database.
        result_func: A function to be applied to the response. Must return an iterable.
        """
        ids, results = self._init_pb(ids, batch_size)

        executor = ThreadPoolExecutor(max_workers=self.num_threads)

        fs = []
        for start in range(0, len(ids), self.batch_size):
            num = len(ids) - start
            num = self.batch_size if num > self.batch_size else num
            f = executor.submit(
                self._elink_batch,
                db=db,
                dbfrom=dbfrom,
                ids=ids[start : start + num],
                result_collector=results,
                result_func=result_func,
                **kwargs,
            )
            fs.append(f)

        futures.wait(fs)
        if self.pbar:
            results.pbar.close()
        return results.results, results.failed

    def efetch(
        self,
        db,
        ids,
        retmode,
        retype,
        batch_size=1,
        result_func=lambda x: [x],
        **kwargs,    ):
        """Interface to the efetch database.
        result_func: A function to be applied to the response. Must return an iterable.
        """
        ids, results = self._init_pb(ids, batch_size)
        executor = ThreadPoolExecutor(max_workers=self.num_threads)

        fs = []
        for start in range(0, len(ids), batch_size):
            num = len(ids) - start
            num = batch_size if num > batch_size else num
            f = executor.submit(
                self._efetch_batch,
                db=db,
                ids=ids[start : start + num],
                result_collector=results,
                result_func=result_func,
                retmode=retmode,
                retype=retype,
                **kwargs,
            )
            fs.append(f)

        futures.wait(fs)
        if self.pbar:
            results.pbar.close()
        return results.results, results.failed

    def esummary(self, db, ids, batch_size=1, result_func=lambda x: [x], **kwargs):
        """Interface to the efetch database.
        result_func: A function to be applied to the response. Must return an iterable.
        """
        ids, results = self._init_pb(ids, batch_size)

        executor = ThreadPoolExecutor(max_workers=self.num_threads)

        fs = []

        for start in range(0, len(ids), batch_size):
            num = len(ids) - start
            num = batch_size if num > batch_size else num
            f = executor.submit(
                self._general_batch,
                db=db,
                ids=ids[start : start + num],
                result_collector=results,
                result_func=result_func,
                emode="esummary",
                **kwargs,
            )
            fs.append(f)

        futures.wait(fs)
        if self.pbar:
            results.pbar.close()
        return results.results, results.failed

    def esearch(self, db, ids, batch_size=1, result_func=lambda x: [x], **kwargs):
        """Interface to the efetch database.
        result_func: A function to be applied to the response. Must return an iterable.
        """

        executor = ThreadPoolExecutor(max_workers=self.num_threads)

        ids, results = self._init_pb(ids, batch_size)

        fs = []
        for start in range(0, len(ids), batch_size):
            num = len(ids) - start
            num = batch_size if num > batch_size else num
            f = executor.submit(
                self._general_batch,
                db=db,
                ids=" OR ".join([str(_) for _ in ids[start : start + num]]),
                result_collector=results,
                result_func=result_func,
                emode="esearch",
                RetMax=self.batch_size * 2,
                **kwargs,
            )
            fs.append(f)

        futures.wait(fs)
        if self.pbar:
            results.pbar.close()
        return results.results, results.failed


if __name__ == "__main__":
    edl = EntrezDownloader(
        # An email address. You might get blocked by the NCBI without specifying one.
        email="l0404th@gmail.com",
        # An API key. You can obtain one by creating an NCBI account. Speeds things up.
        api_key="ccf9847611deebe1446b9814a356f14cde08",
        num_threads=30,  # The number of parallel requests to make
        batch_size=500,  # The number of IDs to fetch per request
        pbar=True,  # Enables a progress bar, requires tqdm package
    )
