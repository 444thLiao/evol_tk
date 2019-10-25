


    results, failed = edl.esearch(db='protein',
                                  ids=remained_records_ids,
                                  result_func=lambda x: Entrez.read(io.StringIO(x))['IdList'])
    all_GI = list(set(results[::]))
    tqdm.write('get pid summary from each one')
    results, failed = edl.esummary(db='protein',
                                   ids=all_GI,
                                   result_func=lambda x: Entrez.read(
                                       io.StringIO(x)))
    id2complete_tax = {}
    for r in tqdm(results):
        aid = r['AccessionVersion'].strip()
        tid = r['TaxId'].real
        lineage = ncbi.get_lineage(tid)
        rank = ncbi.get_rank(lineage)
        rank = {v: k for k, v in rank.items()}
        names = ncbi.get_taxid_translator(lineage)
        rank2names = {k:names.get(tax) for k,tax in rank.items() if k !='no rank'}
        id2complete_tax[aid] = rank2names
    
    b = []
    for _ in SeqIO.parse('./nr_retrieve_removeENV_amoA/filtered_by_kegg.faa',format='fasta'):
        if _.id in id2complete_tax:
            if id2complete_tax[_.id].get('superkingdom','')=='Bacteria' and id2complete_tax[_.id].get('phylum','ENV')!='ENV':
                b.append(_)
    [_ for _ in a if id2complete_tax[_.id].get('superkingdom','')=='Bacteria' and id2complete_tax[_.id].get('phylum','ENV')!='ENV' ]
    # archaea
    len([k for k,v in id2complete_tax.items() if v.get('superkingdom','')=='Archaea'])
    id2tax = {}
    id2org = {}
    for r in tqdm(results):
        aid = r['AccessionVersion']
        tid = r['TaxId'].real
        lineage = ncbi.get_lineage(tid)
        rank = ncbi.get_rank(lineage)
        rank = {v: k for k, v in rank.items()}
        names = ncbi.get_taxid_translator(lineage)
        if names.get(rank.get('phylum',''),'ENV') == 'Proteobacteria':
            id2tax[aid] = names.get(rank.get('class',''),'ENV')
        else:
            id2tax[aid] = names.get(rank.get('phylum',''),'ENV')
        id2org[aid] = names[tid]
    id2tax = {k:v for k,v in id2tax.items() if k in final_ids}
    id2org = {k:v for k,v in id2org.items() if k in final_ids}  