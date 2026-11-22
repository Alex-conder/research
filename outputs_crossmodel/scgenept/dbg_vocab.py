import sys
import numpy as np
import pandas as pd
import anndata as ad

sys.path.insert(0, '/research/scgenept_repo/repo')
from utils.data_loading import match_genes_to_scgpt_vocab_from_adata  # noqa: E402
from utils.scgpt_config import SPECIAL_TOKENS  # noqa: E402

genes = ['POP7', 'AARS1', 'XYZ']
stub = ad.AnnData(X=np.zeros((1, 3), dtype=np.float32),
                  var=pd.DataFrame({'gene_name': genes}))
vocab, gids, dg, g2i = match_genes_to_scgpt_vocab_from_adata(
    '/research/scgenept_repo/repo/models/pretrained/scgpt/vocab.json',
    stub, SPECIAL_TOKENS)
print('dg=', [repr(x) for x in dg])
print('len g2i=', len(g2i))
print('POP7 in g2i:', 'POP7' in g2i)
print('keys sample:', list(g2i.keys())[:3])
print('gids=', gids)
