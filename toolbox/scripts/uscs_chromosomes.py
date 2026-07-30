# This is a dictionary for translating/converting NCBI or ensembl chromosome names to USCS standard
# Central registry mapping species names to their chromosome conversion dicts
CHROM_DICTS = {
    "saccharomyces_cerevisiae": {
        "I": "chrI", "II": "chrII", "III": "chrIII", "IV": "chrIV", 
        "V": "chrV", "VI": "chrVI", "VII": "chrVII", "VIII": "chrVIII", 
        "IX": "chrIX", "X": "chrX", "XI": "chrXI", "XII": "chrXII", 
        "XIII": "chrXIII", "XIV": "chrXIV", "XV": "chrXV", "XVI": "chrXVI", 
        "Mito": "chrM", "MT": "chrM", "chrm": "chrM",
        # NCBI RefSeq Accessions mapping
        "NC_001133.9": "chrI", "NC_001134.8": "chrII", "NC_001135.5": "chrIII",
        "NC_001136.10": "chrIV", "NC_001137.3": "chrV", "NC_001138.5": "chrVI",
        "NC_001139.9": "chrVII", "NC_001140.6": "chrVIII", "NC_001141.2": "chrIX",
        "NC_001142.9": "chrX", "NC_001143.9": "chrXI", "NC_001144.5": "chrXII",
        "NC_001145.3": "chrXIII", "NC_001146.2": "chrXIV", "NC_001147.6": "chrXV",
        "NC_001148.4": "chrXVI", "NC_001224.1": "chrM"
    },

    "caenorhabditis_elegans": {
        # Standard Roman Numerals & Ensembl/WormBase headers
        "I": "chrI", "II": "chrII", "III": "chrIII", 
        "IV": "chrIV", "V": "chrV", "X": "chrX", 
        "MtDNA": "chrM", "MT": "chrM", "Mito": "chrM", "chrm": "chrM",
        # NCBI RefSeq Accessions mapping (WBcel235 / PRJNA13758)
        "NC_003279.8": "chrI", "NC_003280.10": "chrII", "NC_003281.10": "chrIII",
        "NC_003282.8": "chrIV", "NC_003283.11": "chrV", "NC_003284.9": "chrX",
        "NC_001328.1": "chrM"
    },
    
    "homo_sapiens": {
        # Future-proofing: easily populate human chromosome/accession mappings here
        "1": "chr1", "2": "chr2", "MT": "chrM"
    }
}

def get_ucsc_mapping(species_dict):
    """Safely retrieves the dictionary for a species; returns empty dict if missing."""
    return CHROM_DICTS.get(species_dict.lower(), {})