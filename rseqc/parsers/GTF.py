#__version__ = "0.0.2"

import sys
import re

#TODO include support for gziped files

class ParseGTF(object):
    def __init__(self, gene_models, chr = "all"):
        '''This is constructor of ParseGTF that make special gtf object'''

        assert gene_models.split('.')[-1] == 'gtf'

        tags_regex = "\\s.([A-z0-9_.-]+)"
        tags_dict = {
                'gid': 'gene_id',
                'gname': 'gene_name',
                'gtype': 'gene_biotype',
                'tid': 'transcript_id',
                'ttype': 'transcript_biotype',
                'eid': 'exon_id',
                'pid': 'protein_id',
                'tag': 'tag',
                'tsl': 'transcript_support_level'
                }
        
        for key in tags_dict:
                new_value = tags_dict[key]+tags_regex
                tags_dict[key] = new_value

        self.models_dict = {}
        exons_dict = {}
        cds_dict = {}

        go = False
        if chr == "all":
            go = True

        cds_cnter = 0
        pre_gid = None

        with open(gene_models, 'r') as gm:

            for l in gm:
                line = l.strip()
                chk_chr = re.search('(^[A-z0-9_.-]+)\t', line)

                if not line.startswith("#"):
                    if go or chk_chr.group(1) == str(chr):
                        items = line.split('\t')

                        assert len(items) == 9

                        chk_gid = re.search(tags_dict['gid'], items[8])
                        chk_eid = re.search(tags_dict['eid'], items[8])
                        chk_pid = re.search(tags_dict['pid'], items[8])
                        chk_tid = re.search(tags_dict['tid'], items[8])
                        chk_gtype = re.search(tags_dict['gtype'], items[8])
                        chk_ttype = re.search(tags_dict['ttype'], items[8])
                        chk_gname = re.search(tags_dict['gname'], items[8])
                        chk_tag = re.search(tags_dict['tag'], items[8])
                        chk_tsl = re.search(tags_dict['tsl'], items[8])
            
                        if chk_gid.group(1) not in self.models_dict:
                                self.models_dict[chk_gid.group(1)] = {}
                        # in gtf every line must have gene_id
                        gid = chk_gid.group(1)
                        # not checking if mainInfo since if checkExonId isn't null then this is must exist
                        if 'chr' not in self.models_dict[gid]:
                            self.models_dict[gid]['chr'] = items[0]
                        if 'strand' not in self.models_dict[gid]:
                            self.models_dict[gid]['strand'] = items[6]

                        if chk_gtype:
                            if 'gtype' not in self.models_dict[gid]:
                                self.models_dict[gid]['gtype'] = chk_gtype.group(1)
                        else:
                            self.models_dict[gid]['gtype'] = 'NA'
                        # check whether gene name was found and insert it in
                        if chk_gname:
                            if 'gname' not in self.models_dict[gid]:
                                self.models_dict[gid]['gname'] = chk_gname.group(1)
                        else:
                            self.models_dict[gid]['gname'] = 'NA'

                        if chk_tag:
                            if 'tag' not in self.models_dict[gid]:
                                self.models_dict[gid]['tag'] = chk_tag.group(1)
                        else:
                            self.models_dict[gid]['tag'] = 'NA'

                        if items[2].lower() == 'exon':
                            tid = chk_tid.group(1)
                            eid = chk_eid.group(1)
                            # Gathering Transcript and Exons information
                            if 'exons' not in self.models_dict[gid]:
                                self.models_dict[gid]['exons'] = {}
                            self.models_dict[gid]['exons'][eid] = {'start': int(items[3]), 'end': int(items[4])}

                            if gid not in exons_dict:
                                exons_dict[gid] = {}
                            if tid not in exons_dict[gid]:
                                exons_dict[gid][tid] = []
                            exons_dict[gid][tid].append(eid)

                        if items[2].lower() == 'transcript':
                            tid = chk_tid.group(1)

                            if 'tscripts' not in self.models_dict[gid]:
                                self.models_dict[gid]['tscripts'] = {}
                            if tid not in self.models_dict[gid]['tscripts']:
                                self.models_dict[gid]['tscripts'][tid] = {'exons_id': [], 'cds_id': []}

                            attrs = ['ATG', 'TGA', 'tprime_utr', 'fprime_utr']
                            for attr in attrs:
                                if attr not in self.models_dict[gid]['tscripts'][tid]:
                                    self.models_dict[gid]['tscripts'][tid][attr] = {}

                            # make some place holder for certain attributes 
                            if 'tsl' not in self.models_dict[gid]['tscripts'][tid]:
                                self.models_dict[gid]['tscripts'][tid]['tsl'] = 'NA'
                            if 'ttype' not in self.models_dict[gid]['tscripts'][tid]:
                                self.models_dict[gid]['tscripts'][tid]['ttype'] = 'NA'

                            if chk_ttype:
                                self.models_dict[gid]['tscripts'][tid]['ttype'] = chk_ttype.group(1)
                            
                            if chk_tsl:
                                self.models_dict[gid]['tscripts'][tid]['tsl'] = chk_tsl.group(1)

                        if items[2].lower() == 'cds':

                            if pre_gid != gid:
                                cds_cnter = 0
                            cds_cnter += 1
                            pre_gid = gid

                            tid = chk_tid.group(1)
                            pid = chk_pid.group(1)

                            if 'cdss' not in self.models_dict[gid]:
                                self.models_dict[gid]['cdss'] = {}
                            self.models_dict[gid]['cdss'][pid+'.'+str(cds_cnter)] = {'start': int(items[3]), 'end': int(items[4])}

                            if gid not in cds_dict:
                                cds_dict[gid] = {}
                            if tid not in cds_dict[gid]:
                                cds_dict[gid][tid] = []
                            cds_dict[gid][tid].append(pid+'.'+str(cds_cnter))

                            #self.models_dict[gid]['tscripts'][tid]['cds'] = {'start': int(items[3]), 'end': int(items[4])}
                        
                        if items[2].lower() == 'three_prime_utr':
                            self.models_dict[gid]['tscripts'][tid]['tprime_utr'] = {'start': int(items[3]), 'end': int(items[4])}

                        if items[2].lower() == 'five_prime_utr':
                            self.models_dict[gid]['tscripts'][tid]['fprime_utr'] = {'start': int(items[3]), 'end': int(items[4])}

                        if items[2].lower() == 'start_codon':
                            self.models_dict[gid]['tscripts'][tid]['ATG'] = {'start': int(items[3]), 'end': int(items[4])}

                        if items[2].lower() == 'stop_codon':
                            self.models_dict[gid]['tscripts'][tid]['TGA'] = {'start': int(items[3]), 'end': int(items[4])}

        for g in list(exons_dict.keys()):
            for t, v in list(exons_dict[g].items()):
                self.models_dict[g]['tscripts'][t]['exons_id'] = v

        for g in list(cds_dict.keys()):
            for t, v in list(cds_dict[g].items()):
                self.models_dict[g]['tscripts'][t]['cds_id'] = v

class GeneModels(ParseGTF):

    def __init__(self, gene_models, chr = "all"):
        super(GeneModels, self).__init__(gene_models, chr = chr)

        self.genes = []

        for gene, attr in list(self.models_dict.items()):

            edict = attr['exons']
            strand = attr['strand']
            chr = attr['chr']
            gstart = min([v['start'] for v in list(edict.values())])
            gend = max([v['end'] for v in list(edict.values())])

            self.genes.append({'chr': chr,
                               'start': gstart,
                               'end': gend,
                               'id': gene,
                               'strand': strand,
                               'name': attr['gname'],
                               'biotype': attr['gtype']
                               })

    def get_genes(self, chr = 'all', strand = 'both', biotype = 'protein_coding'):
        '''
        Returns a list of lists, where nested list contains seven elements, where each nested list contains gene's elements
       
        Return:

        [.., {'chr': CHR, 'start': START, 'end': END, 'id': ID, 'strand': STRAND, 'name': NAME, 'biotype': BIOTYPE], ..]

        Usage:

        - chr['all'], user can optinally pass in either a string or a list to filter chromosomes of interest

            - chr = 1 or '1' either will work
            - chr = [1, 13, 20] or ['1', '13', '20'] either will work
            - chr = 'all' is a special keyword that returns all chromosomes

        - strand['both'], user can optionally pass in either of the strands in to filter based on gene strand

            - strand = '+'
            - strand = '-'
            - strand = 'both' is a special keyword to get both strands

        - biotype['protein_coding'], user can optionaly pass in either a string or a list 
 
            - biotype = 'snRNA'
            - biotype = ['snRNA', 'lnRNA']
            - biotype = 'all' is a special keyword that returns all biotypes
        '''
        assert isinstance(strand, str)

        tmp = chr
        if chr != 'all' and ( isinstance(chr, str) or isinstance(chr, int) ):
            tmp = [chr] 
        # make sure all entries are converted to strings
        chrs = list(map(str, tmp))

        biotypes = biotype
        if isinstance(biotype, str):
            biotypes = [biotype] 

        if chr == 'all':
            if strand == 'both':
                if biotype == 'all':
                    return self.genes
                else:
                    return [g for g in self.genes if g['biotype'] in biotypes]
            else:
                if biotype == 'all':
                    return [g for g in self.genes if g['strand'] == strand]
                else:
                    return [g for g in self.genes if g['strand'] == strand and g['biotype'] in biotypes]
        else:
            if strand == 'both':
                if biotype == 'all':
                    return [g for g in self.genes if g['chr'] in chrs]
                else:
                    return [g for g in self.genes if g['chr'] in chrs and g['biotype'] in biotypes]
            else:
                if biotype == 'all':
                    return [g for g in self.genes if g['strand'] == strand and g['chr'] in chrs]
                else:
                    return [g for g in self.genes if g['strand'] == strand and g['chr'] in chrs and g['biotype'] in biotypes]

    def get_cds(self):
        '''
        Returns a list of dictionaries, where nested dictionary contains seven elements containing CDS per transcript

        Returns:

        [.., {'chr': CHR, 'start': START, 'end': END, 'id': ID, 'strand': STRAND, 'name': NAME, 'biotype': BIOTYPE}, ..]

        Usage:

        '''

        cds = []

        for gene_obj in self.genes:
            gid = gene_obj['id']
            gname = gene_obj['name']
            chr = gene_obj['chr']
            strand = gene_obj['strand']
            obj = self.models_dict[gid]
            tdict = obj['tscripts']

            for t, v in list(tdict.items()):

                ttype = v['ttype']
                pids = v['cds_id']
                eids = v['exons_id']

                if pids:
                    if strand == '+':
                        pstart = min([obj['cdss'][p]['start'] for p in pids])
                        pend = max([obj['cdss'][p]['end'] for p in pids])
                        cds.append({'chr': chr, 'start': pstart, 'end': pend, 'id': t, 'strand': strand, 'name': gname, 'biotype': ttype})
                    else:
                        pstart = max([obj['cdss'][p]['end'] for p in pids])
                        pend = min([obj['cdss'][p]['start'] for p in pids])
                        cds.append({'chr': chr, 'start': pend, 'end': pstart, 'id': t, 'strand': strand, 'name': gname, 'biotype': ttype})
                else:
                    atg = v['ATG']
                    tga = v['TGA']

                    if atg and tga:
                        if strand == '+':
                            cds.append({'chr': chr, 'start': atg['start'], 'end': tga['end'], 'id': t, 'strand': strand, 'name': gname, 'biotype': ttype})
                        else:
                            cds.append({'chr': chr, 'start': tga['start'], 'end': atg['end'], 'id': t, 'strand': strand, 'name': gname, 'biotype': ttype})
                    elif atg or tga:
                        if atg:
                            if strand == '+':
                                end = max([obj['exons'][e]['end'] for e in eids])
                                cds.append({'chr': chr, 'start': atg['start'], 'end': end, 'id': t, 'strand': strand, 'name': gname, 'biotype': ttype})
                            else:
                                end = min([obj['exons'][e]['start'] for e in eids])
                                cds.append({'chr': chr, 'start': end, 'end': atg['end'], 'id': k, 'strand': strand, 'name': gname, 'biotype': ttype})
                        else:
                            if strand == '+':
                                start = min([obj['exons'][e]['start'] for e in eids])
                                cds.append({'chr': chr, 'start': start, 'end': tga['end'], 'id': t, 'strand': strand, 'name': gname, 'biotype': ttype})
                            else:
                                start = max([obj['exons'][e]['end'] for e in eids])
                                cds.append({'chr': chr, 'start': tga['start'], 'end': start, 'id': k, 'strand': strand, 'name': gname, 'biotype': ttype})
                    else:
                        if strand == '+':
                            start = min([obj['exons'][e]['start'] for e in eids])
                            end = max([obj['exons'][e]['end'] for e in eids])
                            cds.append({'chr': chr, 'start': start, 'end': end, 'id': t, 'strand': strand, 'name': gname, 'biotype': ttype})
                        else:
                            start = max([obj['exons'][e]['end'] for e in eids])
                            end = min([obj['exons'][e]['start'] for e in eids])
                            cds.append({'chr': chr, 'start': end, 'end': start, 'id': t, 'strand': strand, 'name': gname, 'biotype': ttype})

        return cds

    def get_tscripts(self):
        #TODO doesn't look like this method had been finished
        '''
        Returns a list of dictionaries, where nested dictionary contains seven elements per single transcript

        Returns:

        [.., {'chr': CHR, 'start': START, 'end': END, 'id': ID, 'strand': STRAND, 'name': NAME, 'biotype': BIOTYPE], ..]

        Usage:

        '''

        tscripts = []

        for gene_obj in self.genes:
            gid = gene_obj['id']
            gname = gene_obj['name']
            obj = self.models_dict[gid]
            tdict = obj['tscripts']
            tscripts.append(tdict)    
            #for k, v in tdict.items():
            #    self.tscripts.append({'chr': chr, 'start': v['exons_id'], v['end'], k, strand, gene_name, 'NA'])
            #    tscripts.append(v['exons_id'])
        return tscripts

    def get_exons(self):
        '''
        Extract exon regions from GTF file and outputs a list of lists, where nested lists 
        are individual exons and nested list has the following format:
        Chr, Start, End, Id, Strand, Name, Biotype
        '''
        self.exons = []

        #for gene, attr in self.models_dict.items():
        for gene_obj in self.genes:
            gid = gene_obj['chr']
            gname = gene_obj['name']
            obj = self.models_dict[gid]
            edict = obj['exons']
            strand = obj['strand']
            chr = obj['chr']

            for k, v in list(edict.items()):
                self.exons.append([chr, v['start'], v['end'], k, strand, gname, 'NA'])
        return self.exons

    def get_introns(self, biotype = 'protein_coding'):
        '''
        Returns a list of dictionaries, where nested dictionary contains seven elements per single intron

        Returns:

        [.., {'chr': CHR, 'start': START, 'end': END, 'id': ID, 'strand': STRAND, 'name': NAME, 'biotype': BIOTYPE}, ..]

        Usage:

        - biotype['protein_coding'], user can optionaly pass in either a string or a list 

            - biotype = 'snRNA'
            - biotype = ['snRNA', 'lnRNA']

        '''

        biotypes = biotype
        if isinstance(biotype, str):
            biotypes = [biotype] 

        introns = []

        prefix = 'ENSMUSI'
        intron_cnter = 0
  
        for gene_obj in self.genes:
            gid = gene_obj['id']
            gname = gene_obj['name']
            strand = gene_obj['strand']
            chr = gene_obj['chr']
            obj = self.models_dict[gid]
            tdict = obj['tscripts']

            for t, v in list(tdict.items()):
                ttype = v['ttype']
                eids = v['exons_id']

                if ttype in biotypes:

                    starts = [obj['exons'][e]['start'] for e in eids]
                    starts.sort()
                    ends = [obj['exons'][e]['start'] for e in eids]
                    ends.sort()

                    for i in range(len(ends)-1):
                        introns.append({'chr': chr, 'start': ends[i], 'end': starts[i+1], 'id': prefix+str(intron_cnter), 'strand': strand, 'name': gname, 'biotype': ttype})
                        intron_cnter += 1

        return introns

    def get_utrs(self, utr = None, biotype = "protein_coding", tsl = 1):
        '''
        Returns a list of dictionaries, where nested dictionary contains seven elements per single 3 prime UTR

        Returns:

        [.., {'chr': CHR, 'start': START, 'end': END, 'id': ID, 'strand': STRAND, 'name': NAME, 'biotype': BIOTYPE}, ..]

        Usage:

        - tsl[1], user can optinally pass in either a string or a list (tp://asia.ensembl.org/Help/Glossary?id=492)

           - tsl = 1
           - tsl = [1, 2, 3, 4, 5]

        - biotype['protein_coding'], user can optionaly pass in either a string or a list 

            - biotype = 'snRNA'
            - biotype = ['snRNA', 'lnRNA']
        '''

        assert isinstance(utr, int)
        
        if utr == 3:
            utr = 'tprime_utr'
        if utr == 5:
            utr = 'fprime_utr'
        
        utrs = {}
        # this is enable filtering based on biotype
        biotypes = biotype
        if isinstance(biotype, str):
            biotypes = [biotype] 

        tsls = tsl
        if isinstance(tsl, int):
            tsls = [tsl] 
        
        for gene_obj in self.genes:
            gid = gene_obj['id']
            gname = gene_obj['name']
            strand = gene_obj['strand']
            chr = gene_obj['chr']
            obj = self.models_dict[gid]
            tdict = obj['tscripts']

            for t, v in list(tdict.items()):

                ttype = v['ttype']
                eids = v['exons_id']
                pids = v['cds_id']

                if ttype in biotypes:

                    get_tsl = v['tsl']
                    if get_tsl not in utrs:
                        utrs[get_tsl] = []

                    if v[utr]:
                        start = v[utr]['start']
                        end = v[utr]['end']

                        utrs[get_tsl].append({'chr': chr, 'start': start, 'end': end, 'id': t, 'strand': strand, 'name': gname, 'biotype': ttype})

                    else:
                        if utr == 'tprime_utr':
                            tga = v['TGA']
                            # assume 3'UTR for now
                            if tga:
                                if strand == '+':
                                    end = max([obj['exons'][e]['end'] for e in eids])
                                    if tga['end'] != end:
                                        utrs[get_tsl].append({'chr': chr, 'start': tga['end']+1, 'end': end, 'id': t, 'strand': strand, 'name': gname, 'biotype': ttype})
                                else:
                                    end = min([obj['exons'][e]['start'] for e in eids])
                                    if tga['start'] != end:
                                        utrs[get_tsl].append({'chr': chr, 'start': end, 'end': tga['start']-1, 'id': t, 'strand': strand, 'name': gname, 'biotype': ttype})
                                    #TODO might want to output 'NA' for start and end perhaps enabling people to grab transcript that don't have utr
                                    # utrs[get_tsl].append({'chr': chr, 'start': 'NA', 'end': 'NA', 'id': t, 'strand': strand, 'name': gname, 'biotype': ttype})
                            elif pids:
                                if strand == '+':
                                    start = max([obj['cdss'][p]['end'] for p in pids])
                                    end = max([obj['exons'][e]['end'] for e in eids])
                                    if start != end:
                                        utrs[get_tsl].append({'chr': chr, 'start': start+1, 'end': end, 'id': t, 'strand': strand, 'name': gname, 'biotype': ttype})
                                else:
                                    start = min([obj['cdss'][p]['start'] for p in pids])
                                    end = min([obj['exons'][e]['start'] for e in eids])
                                    if start != end:
                                        utrs[get_tsl].append({'chr': chr, 'start': end, 'end': start-1, 'id': t, 'strand': strand, 'name': gname, 'biotype': ttype})

                        if utr == 'fprime_utr':
                            atg = v['ATG']
                            # assume 3'UTR for now
                            if atg:
                                if strand == '+':
                                    end = min([obj['exons'][e]['start'] for e in eids])
                                    if atg['start'] != end:
                                        utrs[get_tsl].append({'chr': chr, 'start': end, 'end': atg['start']-1, 'id': t, 'strand': strand, 'name': gname, 'biotype': ttype})
                                else:
                                    end = max([obj['exons'][e]['end'] for e in eids])
                                    if atg['end'] != end:
                                        utrs[get_tsl].append({'chr': chr, 'start': atg['end']+1, 'end': end, 'id': t, 'strand': strand, 'name': gname, 'biotype': ttype})
                            elif pids:
                                if strand == '+':
                                    start = min([obj['cdss'][p]['start'] for p in pids])
                                    end = min([obj['exons'][e]['start'] for e in eids])
                                    if start != end:
                                        utrs[get_tsl].append({'chr': chr, 'start': end, 'end': start - 1, 'id': t, 'strand': strand, 'name': gname, 'biotype': ttype})
                                else:
                                    start = max([obj['cdss'][p]['end'] for p in pids])
                                    end = max([obj['exons'][e]['end'] for e in eids])
                                    if start != end:
                                        utrs[get_tsl].append({'chr': chr, 'start': start + 1, 'end': end, 'id': t, 'strand': strand, 'name': gname, 'biotype': ttype})


                       #         if int(tga_end) == int(tend):
                       #             if all_tscripts:
                       #                 #utrs[get_tsl].append([chr, 'NA', 'NA', t, strand, gname, ttype])
                       #         else:
                       #             #utrs[get_tsl].append({chr, tga_end+1, tend, t, strand, gname, ttype})
                       #             utrs[get_tsl].append({'chr': chr, 'start': tga_end+1, 'end': tend, 'id': t, 'strand': strand, 'name': gname, 'biotype': ttype})
                       # else:
                       #     tga_end = v['TGA']['start']
                       #     if tga_end != 'NA':
                       #         tend = min([obj['exons'][e]['end'] for e in eids])

                       #         if int(tga_end) == int(tend):
                       #             if all_tscripts:
                       #                 #utrs[get_tsl].append([chr, 'NA', 'NA', t, strand, gname, ttype])
                       #                 utrs[get_tsl].append({'chr': chr, 'start': 'NA', 'end': 'NA', 'id': t, 'strand': strand, 'name': gname, 'biotype': ttype})
                       #         else:
                       #             #utrs[get_tsl].append([chr, tga_end+1, tend, t, strand, gname, ttype])
                       #             utrs[get_tsl].append({'chr': chr, 'start': tga_end-1, 'end': tend, 'id': t, 'strand': strand, 'name': gname, 'biotype': ttype})

        nested_list = [utrs[str(i)] for i in tsls]
        return [item for sublist in nested_list for item in sublist]

    def get_5prime_utrs(self, biotype = "protein_coding", tsl = 1, all_tscripts = False):
        '''
        Returns a list of dictionaries, where nested dictionary contains seven elements per single 5 prime UTR

        Returns:

        [.., {'chr': CHR, 'start': START, 'end': END, 'id': ID, 'strand': STRAND, 'name': NAME, 'biotype': BIOTYPE}, ..]

        Usage:

        - tsl[1], user can optinally pass in either a string or a list (tp://asia.ensembl.org/Help/Glossary?id=492)

           - tsl = 1
           - tsl = [1, 2, 3, 4, 5]

        - biotype['protein_coding'], user can optionaly pass in either a string or a list 

            - biotype = 'snRNA'
            - biotype = ['snRNA', 'lnRNA']

        - all_tscripts[False], user can change this flag if he wants to see transcripts that didn't have annotated 3'UTRs
        '''

        utrs = {}
        # this is enable filtering based on biotype
        biotypes = biotype
        if isinstance(biotype, str):
            biotypes = [biotype] 

        tsls = tsl
        if isinstance(tsl, int):
            tsls = [tsl] 
        
        for gene_obj in self.genes:
            gid = gene_obj['id']
            gname = gene_obj['name']
            strand = gene_obj['strand']
            chr = gene_obj['chr']
            obj = self.models_dict[gid]
            tdict = obj['tscripts']

            for t, v in list(tdict.items()):
                start = v['fprime_utr']['start']
                end = v['fprime_utr']['end']
                ttype = v['ttype']

                if ttype in biotypes:

                    get_tsl = v['tsl']
                    if get_tsl not in utrs:
                        utrs[get_tsl] = []

                    if start != 'NA':
                        utrs[get_tsl].append({'chr': chr, 'start': start, 'end': end, 'id': t, 'strand': strand, 'name': gname, 'biotype': ttype})
                    else:
                        eids = v['exons_id']
                        if strand == '+':
                            utr_end = v['ATG']['start']
                            if utr_end != 'NA':
                                utr_start = min([obj['exons'][e]['end'] for e in eids])
                                
                                if int(utr_start) == int(utr_end):
                                    if all_tscripts:
                                        utrs[get_tsl].append({'chr': chr, 'start': 'NA', 'end': 'NA', 'id': t, 'strand': strand, 'name': gname, 'biotype': ttype})
                                else:
                                    utrs[get_tsl].append({'chr': chr, 'start': utr_start, 'end': utr_end-1, 'id': t, 'strand': strand, 'name': gname, 'biotype': ttype})
                        else:
                            utr_end = v['ATG']['end']
                            if tga_end != 'NA':
                                utr_start = max([obj['exons'][e]['end'] for e in eids])

                                if int(utr_start) == int(utr_end):
                                    if all_tscripts:
                                        utrs[get_tsl].append({'chr': chr, 'start': 'NA', 'end': 'NA', 'id': t, 'strand': strand, 'name': gname, 'biotype': ttype})
                                else:
                                    utrs[get_tsl].append({'chr': chr, 'start': utr_start, 'end': utr_end+1, 'id': t, 'strand': strand, 'name': gname, 'biotype': ttype})

        nested_list = [utrs[str(i)] for i in tsls]
        return [item for sublist in nested_list for item in sublist]
