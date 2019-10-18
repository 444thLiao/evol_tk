# EntrezDownloader

## Usage

You can use the package with the defaults:

```python
from EntrezDownloader import EntrezDownloader

# Create a new downloader instance
edl = EntrezDownloader()

# 100 IDs
ids = [ f'NT_112{num}.1' for num in range(550, 650) ]

# Execute parallel efetch for the specified IDs
results, failed  = edl.efetch(db = 'nuccore', ids = ids)
```

This will generate a list of XML results returned by each batch invocation. Additionally, a list of IDs that were part of failed batches is returned.

**Note that the result order does not necessarily match the original ID order or content. Unmatched IDs might simply be missing. The `failed` list contains only IDs for which the HTTP request against Entrez failed for technical reasons.**

To get more interpretable results, you can pass a function to `efetch()` that processes the results right after they are fetched. E.g you might want to use the Biopython Entrez XML parser to parse the XML and turn it into individual records:

```python
import io
from Bio import Entrez

results, failed = edl.efetch(
  db = 'nuccore',
  ids = ids,
  result_func = lambda xml_text : [ rec for rec in Entrez.parse(io.StringIO(xml_text)) ]
)
```

The function has to return a list. The default function is `lambda xml_text : [xml_text]`.

The `EntrezDownloader` class can be tweaked in the following ways:

```python
edl = EntrezDownloader(
   email = 'foo@bar.baz',                  # An email address. You might get blocked by the NCBI without specifying one.
   api_key = 'abcdefghijklmnopqrstuvwxyz', # An API key. You can obtain one by creating an NCBI account. Speeds things up.
   num_threads = 30,                       # The number of parallel requests to make
   batch_size = 50,                        # The number of IDs to fetch per request
   pbar = True                             # Enables a progress bar, requires tqdm package
   )
```
