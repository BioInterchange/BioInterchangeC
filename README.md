## Table of Contents

1. [Usage](#usage)
   1. [GFF3, GVF, VCF to JSON](#gff3-gvf-vcf-to-json)
   2. [JSON to GFF3, GVF, VCF](#json-to-gff3-gvf-vcf)
   3. [Python API](#python-api)
   4. [MongoDB](#mongodb)
   5. [RethinkDB](#rethinkdb) (legacy support)
2. [Data Model](#data-model)
   1. [Overview](#overview)
      1. [Context Objects](#context-objects)
      2. [Meta Objects](#meta-objects)
      3. [Feature Objects](#feature-objects)
      4. [Summary Objects](#summary-objects)
3. [Python API Reference Cards](#python-api-reference-cards)
4. [JSON Reference Cards](#json-reference-cards)

## Usage

BioInterchange is a command line tool. On OS X, it can be run via the [Terminal](https://en.wikipedia.org/wiki/Terminal_(OS_X)) application. On Linux, well, you use Linux, so you know what a command line is.

Running BioInterchange will cause the software to perform a quick system check before anything else happens. If there are incompatibility problems, then these will be reported and the software exits.

**Note:** The tool BioInterchange 2.0 and its source code are licensed under the short, simple, and permissive [MIT License](https://github.com/indiedotkim/BioInterchangeC/blob/master/LICENSE.txt). If your version of BioInterchange compains about a missing license file, then please update to version 2.0.5 or later.

### Abbreviations

* [GFF3](http://www.sequenceontology.org/gff3.shtml): Generic Feature Format Version 3
* [GVF](http://www.sequenceontology.org/resources/gvf.html): Genome Variation Format
* [JSON](http://json.org/): JavaScript Object Notation
* [JSON-LD](http://json-ld.org/): JSON Linked Data
* [LDJ/LDJSON](https://en.wikipedia.org/wiki/Line_Delimited_JSON): Line Delimited JSON
* [VCF](https://samtools.github.io/hts-specs/VCFv4.2.pdf): Variant Call Format

### GFF3, GVF, VCF to JSON (Converting to a unified data representation.)

Converting genomic file formats to a series of JSON objects:

    biointerchange -o example.ldj example.vcf

Here, the genomics file "example.vcf" is the input. The output of BioInterchange will be written to the file "example.ldj". The "-o" parameter can also be omitted, in which case the output is written directly on the console:

    biointerchange example.vcf

The optional "-u" parameter can be used to add a custom user annotation to the [context object](#context-objects). This can be useful when including BioInterchange in a genomics analysis pipeline:

    biointerchange -u "CNV Analysis; Smith Lab" -o cnv.ldj cnv.gvf

BioInterchange prints its semantic version number with the "-v" parameter:

    biointerchange -v

To show the EULA that came with the software, the "-e" parameter comes in handy:

    biointerchange -e

A brief help text is shown with the "-h" parameter:

    biointerchange -h

### JSON to GFF3, GVF, VCF (Converting to tool-specific data representations.)

Converting JSON objects back to their original genomic file format:

    biointerchange -o example.vcf example.ldj

**Note:** Modified JSON objects (for example, when using the Python API) will only translate back correctly when non-standard key/value pairs are put under the &ldquo;user-defined&rdquo; key.

### Python API

Genomics data can be accessed and processed via the Python API directly through BioInterchange. Each JSON object will be passed on to a Python function and changes made by the Python code will be preserved.

**For example:**

One of BioInterchange's unit-tests calculates the accumulated lengths of all genomic features in a given file. The Python file `simple.py` (see below) in the module `simplepy` can be called like this:

    PYTHONPATH="`pwd`/test-data" ./biointerchange -p simplepy.simple examples/chromosome_BF.gff

This assumes that you have a correct PYTHONHOME set as well. For using BioInterchange as part of a bigger project, it is recommended to add the necessary paths in PYTHONPATH to your shell initialization (e.g., .bashrc, .zshrc).

Just to look at the last LD-JSON line, use this (your JSON pretty-printer might be called something else than `json_pp`):

    PYTHONPATH="`pwd`/test-data" ./biointerchange -p simplepy.simple examples/chromosome_BF.gff | tail -n 1 | json_pp

Accumulated output (see `accumulated-length`):

    {
       "@context" : "https://www.codamono.com/jsonld/biointerchange-s1.json",
       "statistics" : {
          "features" : 432,
          "features-filtered" : 430,
          "meta-lines" : 3,
          "meta-lines-filtered" : false,
          "comment-lines" : 1
       },
       "@type" : "https://www.codamono.com/gfvo-squared#Summary",
       "runtime" : {
          "invocation" : "Wed Mar 10 15:01:28 2021",
          "finish" : "Wed Mar 10 15:01:28 2021",
          "lapsed-seconds" : "0"
       },
       "my-user-data" : {
          "accumulated-length" : 115528
       }
    }

**Compatibility:** Requires Python 3.9.1 or newer. For older versions of Python you might have to recompile BioInterchange from source.

#### Python simple.py

**Directory structure:**

* simplepy/\_\_init\_\_.py
* simplepy/simple.py

**Environment variables:**

* PYTHONHOME needs to be set to where Python is installed on your system
* PYTHONPATH needs to include "./test-data" for this example, so that "simple.py" can be found; the code is in `test-data` because the code is part of the Google Test unit testing for BioInterchange.

**&ldquo;simple.py&rdquo; source code:**

~~~ python
# This variable will be used to accumulate the
# length (in basepairs) of all features that are
# not filtered out in `process_feature` below.
#
# Use this or a similar approach to share data
# across features.
__accumulated_length__ = 0

def setup(context, meta):
    # The context is always being output, but meta
    # could be omitted. How? It is quite simple:
    # returning (a possibly modified) meta Dict
    # will be output, but returning None will
    # output no meta information.
    return meta

def cleanup(summary):
    # Add the accumulated length of all features
    # to the summary:
    summary['my-user-data'] = { 'accumulated-length' : __accumulated_length__ }

    # The Dict `summary` needs to be returned in order
    # to appear in the output. If no summary should appear
    # in the output, then None needs to be returned.
    return summary

def process_feature(feature):
    global __accumulated_length__

    locus = feature['locus']

    length = abs(locus['end'] - locus['start']) + 1

    # If we want to filter out some features (remove them
    # from the output completely), then we can do so my
    # returning None.
    if length < 10000:
        return None

    __accumulated_length__ += length

    feature['my-calculated-length'] = length

    if 'comment' in feature:
        del feature['comment']

    return feature
~~~

### MongoDB

Download and convert cat features:

    wget ftp://ftp.ensembl.org/pub/release-81/variation/vcf/felis_catus/Felis_catus_incl_consequences.vcf.gz
    gunzip Felis_catus_incl_consequences.vcf.gz
    biointerchange -o Felis_catus_incl_consequences.ldj Felis_catus_incl_consequences.vcf

Import line-delimited JSON-LD documents into MongoDB:

    mongoimport --db genomics --collection felis_catus --type json --file Felis_catus_incl_consequences.ldj

Check that the data is actually in the database:

    mongo genomics
    > db.felis_catus.find()

### RethinkDB (legacy support)

**Note:** RethinkDB did not make it as a company, but their database is still available as an open-source project. Adopting it now would probably a bad choice, but if you work with legacy code, then this is how you load/query features in RethinkDB!

Download and convert cat features:

    wget ftp://ftp.ensembl.org/pub/release-81/variation/vcf/felis_catus/Felis_catus_incl_consequences.vcf.gz
    gunzip Felis_catus_incl_consequences.vcf.gz
    biointerchange -o Felis_catus_incl_consequences.ldj Felis_catus_incl_consequences.vcf

Import line-delimited JSON-LD documents into RethinkDB:

    rethinkdb import -f Felis_catus_incl_consequences.ldj --table genomics.felis_catus

Check that the data is actually in the database using RethinkDB's "Data Explorer" (for local installation at `http://localhost:8080/#dataexplorer`):

    r.db('genomics').table('felis_catus')

## Data Model

### Overview

Four kind of JSON objects are created with each execution of BioInterchange -- in this particular order:

1. a "context object"
2. zero, one, or more "meta objects" (pragma and information lines)
3. a series of "feature objects" (the actual genomics feature-data)
4. a "summary object"

Relationships between objects:

* context and summary objectss stand on their own; they are describing how the data was created, how much data was processed, and how long it took
* meta objects themselves are also independent, but they can be referenced (linked to) by feature objects
* feature objects can reference meta objects as well as other feature objects

**In a nutshell:** meta objects and feature objects contain the genomics data from GFF3, GVF, and VCF files; these objects are referencing each other (they are linked). Context and summary objects stand on their own.

*All* objects contain a "@context" key and a "@type" key, which are called context key and type key in the following.

The context key turns the JSON objects into JSON-LD objects. Never heard of JSON-LD? Just skip to the next part and ignore the "@context" key. This is what makes JSON-LD so great: you can handle JSON-LD objects just like JSON objects.

If you do want to make use of the context key, then feed the JSON-LD objects to a Linked Data tool and it will annotate key/value pairs with type information. JSON-LD objects can also be turned into Triple Store compatible data formats, such as RDF N-Triples, RDF N-Quads, and RDF Turtle. Want to see a real-life example of this magic? Head over to the [JSON-LD Playground](http://json-ld.org/playground/) and copy/paste any JSON-LD object into the "JSON-LD Input" text field: the "N-Quads" tab will instantly show a triples representation of the JSON-LD object that can be loaded into a triple store.

**Note:** All JSON-LD annotations are making use of GFVO<sup>2</sup>. Is it legit? Well, it is based on all the lessons learned when designing the <a href="https://peerj.com/articles/933/">Genomic Feature and Variation Ontology</a> (<a href="https://peerj.com/articles/933/">GFVO</a>). So, there is a team behind GFVO<sup>2</sup> too? Yes, stay tuned&hellip;

#### Context Objects

Context objects tell you something about the environment in which the following JSON objects ("meta objects", "feature objects" and "summary objects") were created.

In essence, a context objects contains information about:

* which version of BioInterchange was used
* the input filename
* the filetype of the input (GFF3, GVF, VCF)
* what the output filename was (unless output went to the console)
* whether the Python API was utilized (if so, name of the Python module)
* additional user-defined parameters

**Example:**

    {
        "@context" : "https://www.codamono.com/jsonld/vcf-c1.json",
        "biointerchange-version" : "2.0.0+36",
        "input-filetype" : "VCF",
        "output-file" : null,
        "input-file" : "test-data/playground-vcf.vcf",
        "python-callback" : null,
        "user-defined" : null
    }

Detailed information about each key/value-pair can be found in the JSON Reference Card section.

#### Meta Objects

One meta object is being created for each genomic data file. The meta object contain data of information/pragma lines.

Meta objects are helpful for determining data provenance, establishing links to ontologies that were used, and providing extra annotations that were not attached to features to reduce data redundancy.

For example, the following example provides a textual descriptions for analytic filters that have been applied in a VCF file. Instead of adding this description to every genomic feature for which the filter applies, the textual description is only given once in this meta object, where it can be looked up via the filter keys ("MinAB", "MinDP", "MinMQ", and "Qual").

**Example:**

    {
        "@context" : "https://www.codamono.com/jsonld/vcf-x1.json",
        "vcf-version" : "4.2",
        "reference" : "ftp://ftp-mouse.sanger.ac.uk/ref/GRCm38_68.fa",
        "contig" : {
            "1" : {
                "length" : 195471971
            },
            "3" : {
                "length" : 918278173
            }
        },
        "filter" : {
            "MinAB" : {
                "description" : "Minimum number of alternate bases (INFO/DP4) [5]"
            },
            "MinDP" : {
                "description" : "Minimum read depth (INFO/DP or INFO/DP4) [5]"
            },
            "MinMQ" : {
                "description" : "Minimum RMS mapping quality for SNPs (INFO/MQ) [20]"
            },
            "Qual" : {
                "description" : "Minimum value of the QUAL field [10]"
            }
        },
        "user-defined" : {
            "samtoolsVersion" : [
                "0.1.18-r572"
            ]
        }
    }

#### Feature Objects

Sequences, sequence variations, sequence annotations, genotyping samples, etc., are represented by feature objects.

##### Basic Feature Information

Most basic information about features includes an identifier, a genomic locus, provenane ("SGRP" stands for the Saccharomyces Genome Resequencing Project), feature type information ("SNV" -- single nucleotide variant; a Sequence Ontology term), and references to external databases (SGRP, again, and European Molecular Biology Laboratory accession).

**Example:**

    {
        "@context" : "https://www.codamono.com/jsonld/gvf-f1.json",
        "id" : "76",
        "locus" : {
            "landmark" : "Chr1",
            "start" : 675,
            "end" : 675,
            "strand" : "+"
        },
        "source" : "SGRP",
        "type" : "SNV",
        "dbxref" : [
            "SGRP:s01-675",
            "EMBL:AA816246"
        ]
    }

##### Variant and Reference Sequences

Information about reference sequences is stored under the "reference" key. Variations are stored under the "variants" key and allele specific informations is labeled by "B", "C", etc. ("A" is denoting the reference, which is not explicitly labeled).

**Example:** GVF (basic information omitted)

    {
        "@context" : "https://www.codamono.com/jsonld/gvf-f1.json",
        "reference" : {
            "sequence" : "A",
            "codon" : "GAG"
        },
        "variants" : {
            "B" : {
                "sequence" : "G",
                "codon" : "GAG"
            },
            "C" : {
                "codon" : "GGG",
                "sequence" : "T"
            }
        }
    }

**Example:** VCF (basic information omitted)

    {
        "@context" : "https://www.codamono.com/jsonld/vcf-f1.json",
        "reference" : {
            "sequence" : "G"
        },
        "variants" : {
            "B" : {
                "sequence" : "GA",
                "allele-count" : 36
            }
        }
    }

##### Genomic/Genotyping Samples

VCF genomic files contain information about samples. Reference, variant, and other information is repeated for each sample.

**Example:** VCF (basic and non-sample specific information omitted)

    {
        "@context" : "https://www.codamono.com/jsonld/vcf-f1.json",
        "samples" : [
            {
                "id" : "129P2",
                "variants" : {
                    "B" : {
                        "allele-count-expected" : 4
                    },
                    "C" : {
                        "allele-count-expected" : 5
                    }
                },
                "genotype-quality" : 68,
                "genotype" : {
                    "phased" : false,
                    "alleles" : "BB",
                    "sequences" : [
                        "G",
                        "G"
                    ]
                },
                "AA" : {
                    "genotype-likelihood-phred-scaled" : 53
                },
                "AB" : {
                    "genotype-likelihood-phred-scaled" : 6
                },
                "BB" : {
                    "genotype-likelihood-phred-scaled" : 0
                },
                "user-defined" : {
                    "SP" : 0,
                    "FI" : 0
                }
            }
        ]
    }

#### Summary Objects

Summary objects wrap-up everything that came beforehand. Namely:

* statistics about the genomics data
* runtime information

Genomics data statistics capture how many comments were seen, how much metadata was read, and how much actual genomic features were processed.

Runtime information tell you when BioInterchange was invoked, when it finished processing the data, and how many seconds it took to process the data.

**Example:**

    {
        "@context" : "https://www.codamono.com/jsonld/biointerchange-s1.json",
        "runtime" : {
            "invocation" : "Fri Jul 17 22:34:12 2015",
            "finish" : "Fri Jul 17 22:34:12 2015",
            "lapsed-seconds" : "0"
        },
        "statistics" : {
            "meta-lines" : 19,
            "meta-lines-filtered" : false,
            "features" : 12,
            "features-filtered" : 0,
            "comment-lines" : 0
        }
    }

## Python API Reference Cards

### Function: setup(context, meta)
* **Description:** Initialization function.
* **Parameter &ldquo;context&rdquo;:** Context object as a Python Dict-instance. Modifications to this dict will have no effect on the output of BioInterchange.
* **Parameter &ldquo;meta&rdquo;:** Meta object as a Python Dict-instance.
* **Returns:** The original meta object or a modified version of it. If None is returned, then the output of the meta object will be suppressed by BioInterchange
* **Example:**
<pre>
~~~ python
def setup(context, meta):
    return meta
~~~
</pre>

### Function: cleanup(summary)
* **Description:** Cleanup and finalization function.
* **Parameter &ldquo;summary&rdquo;:** Summary object as a Python Dict-instance.
* **Returns:** The original summary object or a modified version of it. If None is returned, then the output of the summary object will be suppressed by BioInterchange
* **Example:**
<pre>
~~~ python
def cleanup(summary):
    return summary
~~~
</pre>

### Function: process_feature(feature)
* **Description:** Genomic data processing function.
* **Parameter &ldquo;featurey&rdquo;:** Feature object as a Python Dict-instance.
* **Returns:** The original feature object or a modified version of it. If None is returned, then the output of the feature object will be suppressed by BioInterchange
* **Example:**
<pre>
~~~ python
def process_feature(feature):
    return feature
~~~
</pre>

## JSON Reference Cards

### Key: B, C, etc.
* **Description:** Alternative allele information (reference is implicitly labeled "A").
* **Appears in:** feature object
* **JSON Type:** object
* **JSON-LD Type:** gfvo-squared:alleleB, gfvo-squared:alleleC, etc.
* **Genomic Data Source:** GVF (version 1.07), VCF (version 4.2)
* **In Model Version:** 1
* **Example: Variant information on a feature level.**
<pre>
    {
        "@context" : "https://www.codamono.com/jsonld/vcf-f1.json",
        "variants" : {
            "B" : {
                "sequence" : "G"
            }
        }
    }
</pre>
* **Example: Variant information on a sample level.**
<pre>
    {
        "@context" : "https://www.codamono.com/jsonld/vcf-f1.json",
        "samples" : [
            {
                "id" : "129P2",
                "variants" : {
                    "B" : {
                        "allele-count-expected" : 4
                    },
                    "C" : {
                        "allele-count-expected" : 5
                    }
                }
            }
        ]
    }
</pre>

### Key: AA, AB, BB, AC, etc.
* **Description:** Generic genotype information.
* **Appears in:** feature object
* **JSON Type:** object
* **JSON-LD Type:** gfvo-squared:genotypeAA, gfvo-squared:genotypeAB, etc.
* **Genomic Data Source:** GVF (version 1.07), VCF (version 4.2)
* **See Also:** genotype
* **In Model Version:** 1
* **Example:**
<pre>
    {
        "@context" : "https://www.codamono.com/jsonld/vcf-f1.json",
        "samples" : [
            {
                "id" : "129P2",
                "AA" : {
                    "genotype-likelihood-phred-scaled" : 53
                },
                "AB" : {
                    "genotype-likelihood-phred-scaled" : 6
                },
                "BB" : {
                    "genotype-likelihood-phred-scaled" : 0
                },
                "genotype" : {
                    "sequences" : [
                        "G",
                        "A"
                    ],
                    "phased" : false,
                    "alleles" : "AB"
                }
            }
        ]
    }
</pre>

### Key: affected-features
* **Description:** Identifiers of features that are affected by a variant effect.
* **Appears in:** feature object
* **JSON Type:** array of strings
* **JSON-LD Type:** gfvo-squared:affectedFeatures
* **Genomic Data Source:** GVF (version 1.07)
* **In Model Version:** 1
* **Example:**
<pre>
    {
        "@context" : "https://www.codamono.com/jsonld/gvf-f1.json",
        "variants" : {
            "B" : {
                "effects" : [
                    {
                        "affected-features" : [
                            "YAL067W-A"
                        ]
                    }
                ]
            }
        }
    }
</pre>

### Key: affected-feature-type
* **Description:** Type of the features that are affected by a variant effect.
* **Appears in:** feature object
* **JSON Type:** string
* **JSON-LD Type:** gfvo-squared:affectedFeatureType
* **Genomic Data Source:** GVF (version 1.07)
* **In Model Version:** 1
* **Example:**
<pre>
    {
        "@context" : "https://www.codamono.com/jsonld/gvf-f1.json",
        "variants" : {
            "B" : {
                "effects" : [
                    {
                        "affected-feature-type" : "transcript"
                    }
                ]
            }
        }
    }
</pre>

### Key: alias
* **Description:** Aliases of a feature.
* **Appears in:** feature object
* **JSON Type:** array of strings
* **JSON-LD Type:** gfvo-squared:alias
* **Genomic Data Source:** GFF3 (version 1.21), GVF (version 1.07), VCF (version 4.2)
* **Note:** In VCF, only the first entry of the "ID" column becomes an "id" in JSON-LD, whereas second, third, etc., entries are interpreted as "alias" in JSON-LD.
* **In Model Version:** 1
* **Example:**
<pre>
    {
        "@context" : "https://www.codamono.com/jsonld/gff3-f1.json",
        "id" : "rs123",
        "alias" : [
            "feat12",
            "feat12-1"
        ]
    }
</pre>

### Key: alignment
* **Description:** A sequence alignment.
* **Appears in:** feature object
* **JSON Type:** object
* **JSON-LD Type:** gfvo-squared:alignment
* **Genomic Data Source:** GFF3 (version 1.21), GVF (version 1.07), VCF (version 4.2)
* **In Model Version:** 1
* **Example:**
<pre>
    {
        "@context" : "https://www.codamono.com/jsonld/gff-f1.json",
        "alignment" : {
            "id" : "EST23",
            "start" : 1,
            "end" : 21,
            "strand" : null,
            "cigar-string" : "8M3D6M1I6M"
        }
    }
</pre>

### Key: allele-count
* **Description:** Number of alleles in a genotype.
* **Appears in:** feature object
* **JSON Type:** number
* **JSON-LD Type:** gfvo-squared:alleleCount
* **Genomic Data Source:** VCF (version 4.2)
* **In Model Version:** 1
* **Example:**
<pre>
    {
        "@context" : "https://www.codamono.com/jsonld/vcf-f1.json",
        "variants" : {
            "B" : {
                "sequence" : "GA",
                "allele-count" : 36
            }
        }
    }
</pre>

### Key: allele-count-expected
* **Description:** Expected alternate allele count.
* **Appears in:** feature object
* **JSON Type:** number
* **JSON-LD Type:** gfvo-squared:alleleCountExpected
* **Genomic Data Source:** VCF (version 4.2)
* **In Model Version:** 1
* **Example:**
<pre>
    {
        "@context" : "https://www.codamono.com/jsonld/vcf-f1.json",
        "variants" : {
            "B" : {
                "sequence" : "GA",
                "allele-count-expected" : 4
            }
        }
    }
</pre>

### Key: allele-format
* **Description:** Information about key/values used to describe alleles in a data set.
* **Appears in:** meta object
* **JSON Type:** object
* **JSON-LD Type:** gfvo-squared:alleleFormat
* **Genomic Data Source:** VCF (version 4.2)
* **In Model Version:** 1
* **Example:**
<pre>
    {
        "@context" : "https://www.codamono.com/jsonld/vcf-x1.json",
        "allele-format" : {
            "DEL" : {
                "description" : "Deletion"
            }
        }
    }
</pre>

### Key: allele-frequency
* **Description:** Frequency of an allele in a genotype.
* **Appears in:** feature object
* **JSON Type:** number
* **JSON-LD Type:** gfvo-squared:alleleFrequency
* **Genomic Data Source:** VCF (version 4.2)
* **In Model Version:** 1
* **Example:**
<pre>
    {
        "@context" : "https://www.codamono.com/jsonld/vcf-f1.json",
        "variants" : {
            "B" : {
                "sequence" : "GA",
                "allele-frequency" : 1
            }
        }
    }
</pre>

### Key: allele-total-number
* **Description:** Total number of alleles.
* **Appears in:** feature object
* **JSON Type:** number
* **JSON-LD Type:** gfvo-squared:alleleTotalNumber
* **Genomic Data Source:** VCF (version 4.2)
* **In Model Version:** 1
* **Example:**
<pre>
    {
        "@context" : "https://www.codamono.com/jsonld/vcf-f1.json",
        "allele-total-number" : 36
    }
</pre>

### Key: ancestral-allele
* **Description:** Sequence of an ancestral allele.
* **Appears in:** feature object
* **JSON Type:** string
* **JSON-LD Type:** gfvo-squared:ancestralAllele
* **Genomic Data Source:** VCF (version 4.2)
* **In Model Version:** 1
* **Example:**
<pre>
    {
        "@context" : "https://www.codamono.com/jsonld/vcf-f1.json",
        "ancestral-allele" : "C"
    }
</pre>

### Key: annotation-format
* **Description:** Information about keys used to annotate (filter) features in a data set.
* **Appears in:** meta object
* **JSON Type:** object
* **JSON-LD Type:** gfvo-squared:annotationFormat
* **Genomic Data Source:** VCF (version 4.2)
* **In Model Version:** 1
* **Example:**
<pre>
    {
        "@context" : "https://www.codamono.com/jsonld/vcf-x1.json",
        "annotation-format" : {
            "MinDP" : {
                "description" : "Minimum read depth (INFO/DP or INFO/DP4) [5]"
            },
            "Qual" : {
                "description" : "Minimum value of the QUAL field [10]"
            }
        }
    }
</pre>

### Key: annotations
* **Description:** Tags (or labels) assigned to a feature.
* **Appears in:** feature object
* **JSON Type:** array of strings
* **JSON-LD Type:** gfvo-squared:annotations
* **Genomic Data Source:** GVF (version 1.07), VCF (version 4.2)
* **Note:** Annotations are represented by structured pragma statements in GVF and "FILTER" information in VCF.
* **In Model Version:** 1
* **Example: GVF**
<pre>
    {
        "@context" : "https://www.codamono.com/jsonld/gvf-f1.json",
        "annotations" : [
            "SP1"
        ]
    }
</pre>
* **Annotation Source Example (Meta Document):**
<pre>
    {
        "@context" : "https://www.codamono.com/jsonld/gvf-x1.json",
        "technology-platform" : {
            "SP1" : {
                "types" : [
                    "SNV"
                ],
                "comment" : "SNV information documented in Wiki."
            }
        }
    }
</pre>

### Key: attribute-method
* **Description:** Information about attributes in the data set.
* **Appears in:** meta object
* **JSON Type:** object
* **JSON-LD Type:** gfvo-squared:attributeMethod
* **Genomic Data Source:** GVF (version 1.07)
* **In Model Version:** 1
* **Example:**
<pre>
    {
        "@context" : "https://www.codamono.com/jsonld/gvf-x1.json",
        "attribute-method" : {
            "SP7" : {
                "attribute" : "Zygosity",
                "comment" : "Zygosity is reported here as determined in the original study.",
                "sources" : [
                    "SOLiD"
                ],
                "types" : [
                    "SNV"
                ]
            }
        }
    }
</pre>

### Key: average-coverage
* **Description:** Average read coverage.
* **Appears in:** meta object
* **JSON Type:** number
* **JSON-LD Type:** gfvo-squared:averageCoverage
* **Genomic Data Source:** GVF (version 1.07)
* **In Model Version:** 1
* **Example:**
<pre>
    {
        "@context" : "https://www.codamono.com/jsonld/gvf-x1.json",
        "technology-platform" : {
            "SP1" : {
                "average-coverage" : 36
            }
        }
    }
</pre>

### Key: base-quality-rms
* **Description:** Root-mean-square base quality of a genomic position.
* **Appears in:** feature object
* **JSON Type:** number
* **JSON-LD Type:** gfvo-squared:baseQualityRMS
* **Genomic Data Source:** VCF (version 4.2)
* **In Model Version:** 1
* **Example:**
<pre>
    {
        "@context" : "https://www.codamono.com/jsonld/vcf-f1.json",
        "base-quality-rms" : 3
    }
</pre>

### Key: biointerchange-version
* **Description:** <a href="http://semver.org/">Semantic Version</a> number of the BioInterchange software that was used to create the data set.
* **Appears in:** context object
* **JSON Type:** string
* **JSON-LD Type:** gfvo-squared:biointerchangeVersion
* **In Model Version:** 1
* **Example:**
<pre>
    {
        "@context" : "https://www.codamono.com/jsonld/vcf-c1.json",
        "biointerchange-version" : "2.0.0+36"
    }
</pre>

### Key: breakpoint-fasta
* **Description:** A link to a FASTA file that contains sequences specific to breakpoints.
* **Appears in:** meta object
* **JSON Type:** string
* **JSON-LD Type:** gfvo-squared:breakpointFASTA
* **Genomic Data Source:** VCF (version 4.2)
* **In Model Version:** 1
* **Example:**
<pre>
    {
        "@context" : "https://www.codamono.com/jsonld/vcf-x1.json",
        "breakpoint-fasta" : "ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/sv/breakpoint_assemblies.fasta"
    }
</pre>

### Key: build
* **Description:** Name of a genomic build.
* **Appears in:** meta object
* **JSON Type:** string
* **JSON-LD Type:** gfvo-squared:build
* **Genomic Data Source:** GFF3 (version 1.21), GVF (version 1.07)
* **In Model Version:** 1
* **Example:**
<pre>
    {
        "@context" : "https://www.codamono.com/jsonld/gvf-f1.json",
        "genome-build" : [
            {
                "source" : "NCBI",
                "build" : "B36"
            }
        ]
    }
</pre>

### Key: cigar-string
* **Description:** A CIGAR formatted sequence alignment string.
* **Appears in:** feature object
* **JSON Type:** string in CIGAR format
* **JSON-LD Type:** gfvo-squared:cigarString
* **Genomic Data Source:** GFF3 (version 1.21), GVF (version 1.07), VCF (version 4.2)
* **Note:** CIGAR strings from GFF3/GVF are reformatted to match the CIGAR standard: integer followed by a character. When translating JSON-LD back to GFF3/GVF, the alternative GFF3/GVF-specific format is substituted again.
* **In Model Version:** 1
* **Example:**
<pre>
    {
        "@context" : "https://www.codamono.com/jsonld/gff-f1.json",
        "alignment" : {
            "id" : "EST23",
            "start" : 1,
            "end" : 21,
            "strand" : null,
            "cigar-string" : "8M3D6M1I6M"
        }
    }
</pre>

### Key: codon
* **Description:** A codon sequence.
* **Appears in:** feature object
* **JSON Type:** string
* **JSON-LD Type:** gfvo-squared:codon
* **Genomic Data Source:** GVF (version 1.07)
* **In Model Version:** 1
* **Example:**
<pre>
    {
        "@context" : "https://www.codamono.com/jsonld/gvf-f1.json",
        "reference" : {
            "sequence" : "A",
            "codon" : "GAG"
        },
        "variants" : {
            "B" : {
                "sequence" : "G",
                "codon" : "GAG"
            }
        }
    }
</pre>

### Key: codon-phase
* **Description:** Phase of a coding sequence.
* **Appears in:** feature object
* **JSON Type:** number (values: 0, 1, 2)
* **JSON-LD Type:** gfvo-squared:codon-phase
* **Genomic Data Source:** GFF3 (version 1.21)
* **In Model Version:** 1
* **Example:**
<pre>
    {
        "@context" : "https://www.codamono.com/jsonld/gff3-f1.json",
        "codon-phase" : 0
    }
</pre>

### Key: comment
* **Description:** A free-text comment.
* **Appears in:** feature object, meta object
* **JSON Type:** string
* **JSON-LD Type:** gfvo-squared:comment
* **Genomic Data Source:** GFF3 (version 1.21), GVF (version 1.07), VCF (version 4.2)
* **In Model Version:** 1
* **Example:**
<pre>
    {
        "@context" : "https://www.codamono.com/jsonld/gvf-x1.json",
        "technology-platform" : {
            "SP1" : {
                "id" : "SP1",
                "types" : [
                    "SNV"
                ],
                "comment" : "See notes on wiki."
            }
        }
    }
</pre>

### Key: comment-lines
* **Description:** Number of comment lines that were read from a genomics file.
* **Appears in:** summary object
* **JSON Type:** number
* **JSON-LD Type:** gfvo-squared:commentLines
* **In Model Version:** 1
* **Example:**
<pre>
    {
        "@context" : "https://www.codamono.com/jsonld/biointerchange-s1.json",
        "statistics" : {
            "meta-lines" : 19,
            "meta-lines-filtered" : false,
            "features" : 12,
            "features-filtered" : 0,
            "comment-lines" : 0
        }
    }
</pre>

### Key: contig
* **Description:** Information about a continuous sequence region.
* **Appears in:** meta object
* **JSON Type:** object
* **JSON-LD Type:** gfvo-squared:contig
* **Genomic Data Source:** GFF3 (version 1.21), GVF (version 1.07), VCF (version 4.2)
* **In Model Version:** 1
* **Example:**
<pre>
    {
        "@context" : "https://www.codamono.com/jsonld/gff3-x1.json",
        "contig" : {
            "chr3" : {
                "start" : 1,
                "end" : 99,
                "length" : 99
            }
        }
    }
</pre>

### Key: data-source
* **Description:** Information about the source of genomic/proteomic features.
* **Appears in:** meta object
* **JSON Type:** object
* **JSON-LD Type:** gfvo-squared:dataSource
* **Genomic Data Source:** GVF (version 1.07)
* **In Model Version:** 1
* **Example:**
<pre>
    {
        "@context" : "https://www.codamono.com/jsonld/gvf-x1.json",
        "data-source" : {
            "SP3" : {
                "data-type" : "DNA sequence",
                "comment" : "NCBI Short Read Archive (http://www.ncbi.nlm.nih.gov/Traces/sra)",
                "types" : [
                    "SNV"
                ],
                "dbxref" : [
                    "SRA:SRA008175"
                ],
                "sources" : [
                    "MAQ"
                ]
            }
        }
    }
</pre>

### Key: data-type
* **Description:** Type of data that is presented by a data source.
* **Appears in:** meta object
* **JSON Type:** string
* **JSON-LD Type:** gfvo-squared:dataType
* **Genomic Data Source:** GVF (version 1.07)
* **See Also:** Key: data-source
* **In Model Version:** 1
* **Example:**
<pre>
    {
        "@context" : "https://www.codamono.com/jsonld/gvf-x1.json",
        "data-source" : {
            "SP3" : {
                "data-type" : "DNA sequence",
                "comment" : "NCBI Short Read Archive (http://www.ncbi.nlm.nih.gov/Traces/sra)",
                "types" : [
                    "SNV"
                ],
                "dbxref" : [
                    "SRA:SRA008175"
                ],
                "sources" : [
                    "MAQ"
                ]
            }
        }
    }
</pre>

### Key: dbxref
* **Description:** External database references.
* **Appears in:** feature object, meta object
* **JSON Type:** string
* **JSON-LD Type:** gfvo-squared:dbxref
* **Genomic Data Source:** GFF3 (version 1.21), GVF (version 1.07)
* **In Model Version:** 1
* **Example:**
<pre>
    {
        "@context" : "https://www.codamono.com/jsonld/gvf-f1.json",
        "dbxref" : [
            "SGRP:s01-675",
            "EMBL:AA816246"
        ]
    }
</pre>

### Key: depth
* **Description:** Read depth.
* **Appears in:** feature object
* **JSON Type:** number
* **JSON-LD Type:** gfvo-squared:depth
* **Genomic Data Source:** VCF (version 4.2)
* **In Model Version:** 1
* **Example:**
<pre>
    {
        "@context" : "https://www.codamono.com/jsonld/vcf-f1.json",
        "depth" : 75,
        "samples" : [
            {
                "id" : "129P2",
                "depth" : 2
            }
        ]
    }
</pre>

### Key: effect
* **Description:** Name of the effect a variant has on another feature.
* **Appears in:** feature object
* **JSON Type:** string
* **JSON-LD Type:** gfvo-squared:effect
* **Genomic Data Source:** GVF (version 1.07)
* **In Model Version:** 1
* **Example:**
<pre>
    {
        "@context" : "https://www.codamono.com/jsonld/gvf-f1.json",
        "variants" : {
            "B" : {
                "effects" : [
                    {
                        "effect" : "upstream_gene_variant"
                    }
                ]
            }
        }
    }
</pre>

### Key: effects
* **Description:** Container for variant effects.
* **Appears in:** feature object
* **JSON Type:** array of objects
* **JSON-LD Type:** gfvo-squared:effects
* **Genomic Data Source:** GVF (version 1.07)
* **In Model Version:** 1
* **Example:**
<pre>
    {
        "@context" : "https://www.codamono.com/jsonld/gvf-f1.json",
        "variants" : {
            "B" : {
                "effects" : [
                    {
                        "effect" : "upstream_gene_variant"
                    }
                ]
            }
        }
    }
</pre>

### Key: end
* **Description:** End coordinate of a feature or alignment.
* **Appears in:** feature object
* **JSON Type:** number
* **JSON-LD Type:** gfvo-squared:end
* **Genomic Data Source:** GFF3 (version 1.21), GVF (version 1.07), VCF (version 4.2)
* **In Model Version:** 1
* **Example:**
<pre>
    {
        "@context" : "https://www.codamono.com/jsonld/gff3-f1.json",
        "locus" : {
            "landmark" : "Chr1",
            "start" : 290207,
            "end" : 291002,
            "strand" : "+"
        },
        "alignment" : {
            "id" : "EST23",
            "start" : 1,
            "end" : 21,
            "strand" : null,
            "cigar-string" : "8M3D6M1I6M"
        }
    }
</pre>

### Key: experimentally-validated
* **Description:** Indicates whether sequence variant has been experimentally validated.
* **Appears in:** feature object
* **JSON Type:** boolean
* **JSON-LD Type:** gfvo-squared:experimentallyValidated
* **Genomic Data Source:** VCF (version 4.2)
* **In Model Version:** 1
* **Example:**
<pre>
    {
        "@context" : "https://www.codamono.com/jsonld/vcf-f1.json",
        "experimentally-validated" : false
    }
</pre>

### Key: feature-format
* **Description:** Information about key/values used to describe feature-centric data in a data set.
* **Appears in:** meta object
* **JSON Type:** object
* **JSON-LD Type:** gfvo-squared:featureFormat
* **Genomic Data Source:** VCF (version 4.2)
* **In Model Version:** 1
* **Example:**
<pre>
    {
        "@context" : "https://www.codamono.com/jsonld/vcf-x1.json",
        "feature-format" : {
            "membership-hapmap-2" : {
                "type" : "String",
                "number" : 1,
                "description" : "HapMap 2 membership"
            },
            "AF1" : {
                "type" : "Float",
                "number" : 1,
                "description" : "Max-likelihood estimate of the first ALT allele frequency (assuming HWE)"
            }
        }
    }
</pre>

### Key: features
* **Description:** Number of genomic-features lines that were read from a genomics file.
* **Appears in:** summary object
* **JSON Type:** number
* **JSON-LD Type:** gfvo-squared:features
* **In Model Version:** 1
* **Example:**
<pre>
    {
        "@context" : "https://www.codamono.com/jsonld/biointerchange-s1.json",
        "statistics" : {
            "meta-lines" : 19,
            "meta-lines-filtered" : false,
            "features" : 12,
            "features-filtered" : 0,
            "comment-lines" : 0
        }
    }
</pre>

### Key: features-filtered
* **Description:** Number of genomic-features lines that were filtered through the Python API.
* **Appears in:** summary object
* **JSON Type:** number
* **JSON-LD Type:** gfvo-squared:featuresFiltered
* **In Model Version:** 1
* **Example:**
<pre>
    {
        "@context" : "https://www.codamono.com/jsonld/biointerchange-s1.json",
        "statistics" : {
            "meta-lines" : 19,
            "meta-lines-filtered" : false,
            "features" : 12,
            "features-filtered" : 0,
            "comment-lines" : 0
        }
    }
</pre>

### Key: file-date
* **Description:** Creation date of the file that contains the represented genomics data.
* **Appears in:** meta object
* **JSON Type:** string
* **JSON-LD Type:** gfvo-squared:fileDate
* **Genomic Data Source:** GVF (version 1.07), VCF (version 4.2)
* **In Model Version:** 1
* **Example:**
<pre>
    {
        "@context" : "https://www.codamono.com/jsonld/gvf-x1.json",
        "file-date" : "2015-03-08"
    }
</pre>

### Key: finish
* **Description:** Finish time of the BioInterchange software.
* **Appears in:** summary object
* **JSON Type:** string
* **JSON-LD Type:** gfvo-squared:finish
* **In Model Version:** 1
* **Example:**
<pre>
    {
        "@context" : "https://www.codamono.com/jsonld/biointerchange-s1.json",
         "runtime" : {
            "invocation" : "Fri Jul 17 17:00:54 2015",
            "finish" : "Fri Jul 17 17:00:56 2015",
            "lapsed-seconds" : "2"
        }
    }
</pre>

### Key: genome-build
* **Description:** Information about the underlying genome builds of a data set.
* **Appears in:** meta object
* **JSON Type:** array of objects
* **JSON-LD Type:** gfvo-squared:genomeBuild
* **Genomic Data Source:** GFF3 (version 1.21), GVF (version 1.07)
* **In Model Version:** 1
* **Example:**
<pre>
    {
        "@context" : "https://www.codamono.com/jsonld/gvf-f1.json",
        "genome-build" : [
            {
                "source" : "NCBI",
                "build" : "B36"
            }
        ]
    }
</pre>

### Key: genomic-source
* **Description:** Information about the genomic origin source (LOINC code) of feature data.
* **Appears in:** meta object
* **JSON Type:** string
* **JSON-LD Type:** gfvo-squared:genomicSource
* **Genomic Data Source:** GVF (version 1.07)
* **In Model Version:** 1
* **Example:**
<pre>
    {
        "@context" : "https://www.codamono.com/jsonld/gvf-x1.json",
        "genomic-source" : "somatic"
    }
</pre>

### Key: genotype
* **Description:** Feature or sample specific genotype information.
* **Appears in:** feature object
* **JSON Type:** object
* **JSON-LD Type:** gfvo-squared:genotype
* **Genomic Data Source:** GVF (version 1.07), VCF (version 4.2)
* **See Also:** AA, AB, etc.
* **In Model Version:** 1
* **Example:**
<pre>
    {
        "@context" : "https://www.codamono.com/jsonld/vcf-f1.json",
        "samples" : [
            {
                "id" : "129P2",
                "AA" : {
                    "genotype-likelihood-phred-scaled" : 53
                },
                "AB" : {
                    "genotype-likelihood-phred-scaled" : 6
                },
                "BB" : {
                    "genotype-likelihood-phred-scaled" : 0
                },
                "genotype" : {
                    "sequences" : [
                        "G",
                        "A"
                    ],
                    "phased" : false,
                    "alleles" : "AB"
                }
            }
        ]
    }
</pre>

### Key: genotype-format
* **Description:** Information about keys/values used with genotypes in a data set.
* **Appears in:** meta object
* **JSON Type:** object
* **JSON-LD Type:** gfvo-squared:genotypeFormat
* **Genomic Data Source:** VCF (version 4.2)
* **In Model Version:** 1
* **Example:**
<pre>
    {
        "@context" : "https://www.codamono.com/jsonld/vcf-x1.json",
        "genotype-format" : {
            "depth" : {
                "type" : "Integer".
                "number" : 1,
                "description" : "# high-quality bases"
            },
            "genotype-likelihood" : {
                "type" : "Float",
                "number" : 3,
                "description" : "Likelihoods for RR,RA,AA genotypes (R=ref,A=alt)"
            },
            "SP" : {
                "type" : "Integer",
                "number" : 1,
                "description" : "Phred-scaled strand bias P-value"
            }
        }
    }
</pre>

### Key: genotype-likelihood
* **Description:** Log-10 scaled genotype likelihood.
* **Appears in:** feature object
* **JSON Type:** number
* **JSON-LD Type:** gfvo-squared:genotypeLikelihood
* **Genomic Data Source:** VCF (version 4.2)
* **See Also:** AA, AB, etc.
* **In Model Version:** 1
* **Example:**
<pre>
    {
        "@context" : "https://www.codamono.com/jsonld/vcf-f1.json",
        "samples" : [
            {
                "id" : "129P2",
                "AA" : {
                    "genotype-likelihood" : 12
                },
                "AB" : {
                    "genotype-likelihood" : 3
                },
                "BB" : {
                    "genotype-likelihood" : 0
                },
                "genotype" : {
                    "sequences" : [
                        "G",
                        "A"
                    ],
                    "phased" : false,
                    "alleles" : "AB"
                }
            }
        ]
    }
</pre>

### Key: genotype-likelihood-phred-scaled
* **Description:** Phred-scaled genotype likelihood.
* **Appears in:** feature object
* **JSON Type:** number
* **JSON-LD Type:** gfvo-squared:genotypeLikelihoodPhredScaled
* **Genomic Data Source:** VCF (version 4.2)
* **See Also:** AA, AB, etc.
* **In Model Version:** 1
* **Example:**
<pre>
    {
        "@context" : "https://www.codamono.com/jsonld/vcf-f1.json",
        "samples" : [
            {
                "id" : "129P2",
                "AA" : {
                    "genotype-likelihood-phred-scaled" : 53
                },
                "AB" : {
                    "genotype-likelihood-phred-scaled" : 12
                },
                "BB" : {
                    "genotype-likelihood-phred-scaled" : 0
                },
                "genotype" : {
                    "sequences" : [
                        "G",
                        "A"
                    ],
                    "phased" : false,
                    "alleles" : "AB"
                }
            }
        ]
    }
</pre>

### Key: gff-version
* **Description:** Genomics file-format versioning of GFF3 data sources.
* **Appears in:** meta object
* **JSON Type:** string
* **JSON-LD Type:** gfvo-squared:gffVersion
* **Genomic Data Source:** GFF3 (version 1.21), GVF (version 1.07)
* **In Model Version:** 1
* **Example:**
<pre>
    {
        "@context" : "https://www.codamono.com/jsonld/gff3-x1.json",
        "gff-version" : "1.21"
    }
</pre>

### Key: global
* **Description:** Global assignment that applies to all features in the data set.
* **Appears in:** meta object
* **JSON Type:** object
* **JSON-LD Type:** gfvo-squared:global
* **Genomic Data Source:** GFF3 (version 1.21), GVF (version 1.07), VCF (version 4.2)
* **In Model Version:** 1
* **Example:**
<pre>
    {
        "@context" : "https://www.codamono.com/jsonld/gvf-x1.json",
        "technology-platform" : {
            "global" : {
                "comment" : "Preliminary data."
            }
        }
    }
</pre>

### Key: gvf-version
* **Description:** Genomics file-format versioning of GVF data sources.
* **Appears in:** meta object
* **JSON Type:** string
* **JSON-LD Type:** gfvo-squared:gvfVersion
* **Genomic Data Source:** GVF (version 1.07)
* **In Model Version:** 1
* **Example:**
<pre>
    {
        "@context" : "https://www.codamono.com/jsonld/gvf-x1.json",
        "gvf-version" : "1.07"
    }
</pre>

### Key: id
* **Description:** An identifier.
* **Appears in:** feature object, meta object
* **JSON Type:** string
* **JSON-LD Type:** gfvo-squared:id
* **Genomic Data Source:** GFF3 (version 1.21), GVF (version 1.07), VCF (version 4.2)
* **Note:** Multiple "ID" entries in VCF: only the first entry becomes an "id" in JSON-LD, where as the remaining identifiers become "alias" entries; this behaviour is in concordance with GFF3/GVF representations of identifiers.
* **In Model Version:** 1
* **Example:**
<pre>
    {
        "@context" : "https://www.codamono.com/jsonld/gff3-f1.json",
        "id" : "ENSG00000139618"
    }
</pre>

### Key: individuals
* **Description:** Identifiers of sequenced individuals.
* **Appears in:** meta object
* **JSON Type:** array of strings
* **JSON-LD Type:** gfvo-squared:individuals
* **Genomic Data Source:** GVF (version 1.07)
* **In Model Version:** 1
* **Example:**
<pre>
    {
        "@context" : "https://www.codamono.com/jsonld/gvf-x1.json",
        "individuals" : [
            "NA18507",
            "NA12878",
            "NA19240"
        ]
    }
</pre>

### Key: input-file
* **Description:** Name and relative path of the input file from which BioInterchange read data.
* **Appears in:** context object
* **JSON Type:** string
* **JSON-LD Type:** gfvo-squared:inputFile
* **In Model Version:** 1
* **Example:**
<pre>
    {
        "@context" : "https://www.codamono.com/jsonld/vcf-c1.json",
        "input-file" : "test-data/playground-vcf.vcf"
    }
</pre>

### Key: input-filetype
* **Description:** File type of the input.
* **Appears in:** context object
* **JSON Type:** string (either &ldquo;GFF3&rdquo;, &ldquo;GVF&rdquo;, or &ldquo;VCF&rdquo;)
* **JSON-LD Type:** gfvo-squared:inputFiletype
* **In Model Version:** 1
* **Example:**
<pre>
    {
        "@context" : "https://www.codamono.com/jsonld/vcf-c1.json",
        "input-filetype" : "VCF"
    }
</pre>

### Key: invocation
* **Description:** Invocation time of the BioInterchange software.
* **Appears in:** summary object
* **JSON Type:** string
* **JSON-LD Type:** gfvo-squared:invocation
* **In Model Version:** 1
* **Example:**
<pre>
    {
        "@context" : "https://www.codamono.com/jsonld/biointerchange-s1.json",
         "runtime" : {
            "invocation" : "Fri Jul 17 17:00:54 2015",
            "finish" : "Fri Jul 17 17:00:56 2015",
            "lapsed-seconds" : "2"
        }
    }
</pre>

### Key: landmark
* **Description:** A named genomic or proteomic landmark.
* **Appears in:** feature object
* **JSON Type:** string
* **JSON-LD Type:** gfvo-squared:landmark
* **Genomic Data Source:** GFF3 (version 1.21), GVF (version 1.07), VCF (version 4.2)
* **In Model Version:** 1
* **Example:**
<pre>
    {
        "@context" : "https://www.codamono.com/jsonld/gff3-f1.json",
        "locus" : {
            "landmark" : "Chr1",
            "start" : 290207,
            "end" : 291002
        }
    }
</pre>

### Key: landmarks
* **Description:** List of landmark identifiers.
* **Appears in:** meta object
* **JSON Type:** array of strings
* **JSON-LD Type:** gfvo-squared:landmarks
* **Genomic Data Source:** GVF (version 1.07)
* **In Model Version:** 1
* **Example:**
<pre>
    {
        "@context" : "https://www.codamono.com/jsonld/gvf-f1.json",
        "technology-platform" : [
            {
                "id" : "SP1",
                "landmarks" : [
                    "Chr1", "ChrX"
                ]
            }
        ]
    }
</pre>

### Key: lapsed-seconds
* **Description:** Number of seconds that the BioInterchange software took for one execution.
* **Appears in:** summary object
* **JSON Type:** string
* **JSON-LD Type:** gfvo-squared:lapsed-seconds
* **In Model Version:** 1
* **Example:**
<pre>
    {
        "@context" : "https://www.codamono.com/jsonld/biointerchange-s1.json",
         "runtime" : {
            "invocation" : "Fri Jul 17 17:00:54 2015",
            "finish" : "Fri Jul 17 17:00:56 2015",
            "lapsed-seconds" : "2"
        }
    }
</pre>

### Key: length
* **Description:** Length of a genomic feature or continuous sequence region.
* **Appears in:** meta object
* **JSON Type:** number
* **JSON-LD Type:** gfvo-squared:length
* **Genomic Data Source:** GFF3 (version 1.21), GVF (version 1.07), VCF (version 4.2)
* **In Model Version:** 1
* **Example:**
<pre>
    {
        "@context" : "https://www.codamono.com/jsonld/gff3-x1.json",
        "contig" : {
            "chr3" : {
                "start" : 1,
                "end" : 99,
                "length" : 99
            }
        }
    }
</pre>

### Key: locus
* **Description:** A genomic or proteomic locus.
* **Appears in:** feature object
* **JSON Type:** object
* **JSON-LD Type:** gfvo-squared:locus
* **Genomic Data Source:** GFF3 (version 1.21), GVF (version 1.07), VCF (version 4.2)
* **In Model Version:** 1
* **Example:**
<pre>
    {
        "@context" : "https://www.codamono.com/jsonld/gff3-f1.json",
        "locus" : {
            "landmark" : "Chr1",
            "start" : 290207,
            "end" : 291002
        }
    }
</pre>

### Key: mapping-quality-rms
* **Description:** Root-mean-square mapping quality.
* **Appears in:** feature object
* **JSON Type:** number
* **JSON-LD Type:** gfvo-squared:mappingQualityRMS
* **Genomic Data Source:** VCF (version 4.2)
* **In Model Version:** 1
* **Example:**
<pre>
    {
        "@context" : "https://www.codamono.com/jsonld/vcf-f1.json",
        "mapping-quality-rms" : 29
    }
</pre>

### Key: membership-1000G
* **Description:** Membership in the 1000 Genomes (1000G) project.
* **Appears in:** feature object
* **JSON Type:** boolean
* **JSON-LD Type:** gfvo-squared:membership1000G
* **Genomic Data Source:** VCF (version 4.2)
* **In Model Version:** 1
* **Example:**
<pre>
    {
        "@context" : "https://www.codamono.com/jsonld/vcf-f1.json",
        "membership-1000G" : true
    }
</pre>

### Key: membership-dbsnp
* **Description:** Membership in dbSNP.
* **Appears in:** feature object
* **JSON Type:** boolean
* **JSON-LD Type:** gfvo-squared:membershipDbSNP
* **Genomic Data Source:** VCF (version 4.2)
* **In Model Version:** 1
* **Example:**
<pre>
    {
        "@context" : "https://www.codamono.com/jsonld/vcf-f1.json",
        "membership-dbsnp" : true
    }
</pre>

### Key: membership-hapmap-2
* **Description:** Membership in HapMap 2.
* **Appears in:** feature object
* **JSON Type:** boolean
* **JSON-LD Type:** gfvo-squared:membershipHapMap2
* **Genomic Data Source:** VCF (version 4.2)
* **In Model Version:** 1
* **Example:**
<pre>
    {
        "@context" : "https://www.codamono.com/jsonld/vcf-f1.json",
        "membership-hapmap-2" : true
    }
</pre>

### Key: membership-hapmap-3
* **Description:** Membership in HapMap 3.
* **Appears in:** feature object
* **JSON Type:** boolean
* **JSON-LD Type:** gfvo-squared:membershipHapMap3
* **Genomic Data Source:** VCF (version 4.2)
* **In Model Version:** 1
* **Example:**
<pre>
    {
        "@context" : "https://www.codamono.com/jsonld/vcf-f1.json",
        "membership-hapmap-3" : true
    }
</pre>

### Key: meta-lines
* **Description:** Number of genomic-meta/pragma lines that were read from a genomics file.
* **Appears in:** summary object
* **JSON Type:** number
* **JSON-LD Type:** gfvo-squared:metaLines
* **In Model Version:** 1
* **Example:**
<pre>
    {
        "@context" : "https://www.codamono.com/jsonld/biointerchange-s1.json",
        "statistics" : {
            "meta-lines" : 19,
            "meta-lines-filtered" : false,
            "features" : 12,
            "features-filtered" : 0,
            "comment-lines" : 0
        }
    }
</pre>

### Key: meta-lines-filtered
* **Description:** Number of genomic-features lines that were filtered through the Python API.
* **Appears in:** summary object
* **JSON Type:** number
* **JSON-LD Type:** gfvo-squared:metaLinesFiltered
* **In Model Version:** 1
* **Example:**
<pre>
    {
        "@context" : "https://www.codamono.com/jsonld/biointerchange-s1.json",
        "statistics" : {
            "meta-lines" : 19,
            "meta-lines-filtered" : false,
            "features" : 12,
            "features-filtered" : 0,
            "comment-lines" : 0
        }
    }
</pre>

### Key: ontology
* **Description:** A reference to an external ontology.
* **Appears in:** meta object
* **JSON Type:** string
* **JSON-LD Type:** gfvo-squared:ontology
* **Genomic Data Source:** GFF3 (version 1.21), GVF (version 1.07)
* **In Model Version:** 1
* **Example:**
<pre>
    {
        "@context" : "https://www.codamono.com/jsonld/gvf-x1.json",
        "phenotype-description" : {
            "global" : {
                "ontology" : "http://www.human-phenotype-ontology.org/human-phenotype-ontology.obo.gz"
            }
        }
    }
</pre>

### Key: ontology-term
* **Description:** Ontology terms associated with a genomic feature.
* **Appears in:** feature object
* **JSON Type:** array of strings
* **JSON-LD Type:** gfvo-squared:ontology-term
* **Genomic Data Source:** GFF3 (version 1.21), GVF (version 1.07)
* **In Model Version:** 1
* **Example:**
<pre>
    {
        "@context" : "https://www.codamono.com/jsonld/gff3-f1.json",
        "ontology-term" : [
            "GO:0046703"
        ]
    }
</pre>

### Key: output-file
* **Description:** Filename and relative path to the output of the BioInterchange software.
* **Appears in:** context object
* **JSON Type:** string (&ldquo;null&rdquo;, if output was written to the console)
* **JSON-LD Type:** gfvo-squared:outputFile
* **In Model Version:** 1
* **Example:**
<pre>
    {
        "@context" : "https://www.codamono.com/jsonld/vcf-c1.json",
        "output-file" : "example.ldj"
    }
</pre>

### Key: phase-set
* **Description:** Phase set identifier; indicates to which phase set a phased genotype belongs.
* **Appears in:** feature object
* **JSON Type:** number
* **JSON-LD Type:** gfvo-squared:phaseSet
* **Genomic Data Source:** VCF (version 4.2)
* **In Model Version:** 1
* **Example:**
<pre>
    {
        "@context" : "https://www.codamono.com/jsonld/vcf-f1.json",
        "samples" : [
            {
                "id" : "129P2",
                "genotype" : {
                    "alleles" : "BB",
                    "phased" : true,
                    "sequences" : [
                        "A",
                        "A"
                    ]
                },
                "phase-set" : 1
            }
        ]
    }
</pre>

### Key: phased
* **Description:** Indicates whether a genotype is phased or unphased.
* **Appears in:** feature object
* **JSON Type:** boolean
* **JSON-LD Type:** gfvo-squared:phased
* **Genomic Data Source:** GVF (version 1.07), VCF (version 4.2)
* **In Model Version:** 1
* **Example:**
<pre>
    {
        "@context" : "https://www.codamono.com/jsonld/vcf-f1.json",
        "samples" : [
            {
                "id" : "129P2",
                "genotype" : {
                    "alleles" : "BB",
                    "phased" : true,
                    "sequences" : [
                        "A",
                        "A"
                    ]
                },
                "phase-set" : 1
            }
        ]
    }
</pre>

### Key: phased-genotypes
* **Description:** Marks phased genotypes and provides extra information about them.
* **Appears in:** meta object
* **JSON Type:** boolean
* **JSON-LD Type:** gfvo-squared:phasedGenotypes
* **Genomic Data Source:** GVF (version 1.07)
* **In Model Version:** 1
* **Example:**
<pre>
    {
        "@context" : "https://www.codamono.com/jsonld/gvf-x1.json",
        "phased-genotypes" : {
            "SP11" : {
                "types" : [
                    "SNV"
                ]
            }
        }
    }
</pre>

### Key: phasing-quality
* **Description:** Phred-scaled probability that alleles are ordered incorrectly.
* **Appears in:** feature object
* **JSON Type:** number
* **JSON-LD Type:** gfvo-squared:phasingQuality
* **Genomic Data Source:** VCF (version 4.2)
* **In Model Version:** 1
* **Example:**
<pre>
    {
        "@context" : "https://www.codamono.com/jsonld/vcf-f1.json",
        "samples" : [
            {
                "id" : "129P3",
                "genotype" : {
                    "alleles" : "AB",
                    "phased" : true,
                    "sequences" : [
                        "A",
                        "T"
                    ]
                },
                "phasing-quality" : 4
            }
        ]
    }
</pre>

### Key: phenotype-description
* **Description:** Additional information about phenotypes.
* **Appears in:** meta object
* **JSON Type:** object
* **JSON-LD Type:** gfvo-squared:phenotypeDescription
* **Genomic Data Source:** GVF (version 1.07)
* **In Model Version:** 1
* **Example:**
<pre>
    {
        "@context" : "https://www.codamono.com/jsonld/gvf-x1.json",
        "phenotype-description" : {
            "global" : {
                "ontology" : "http://www.human-phenotype-ontology.org/human-phenotype-ontology.obo.gz"
            }
        }
    }
</pre>

### Key: population
* **Description:** Population code that is assigned to an individual.
* **Appears in:** meta object
* **JSON Type:** string
* **JSON-LD Type:** gfvo-squared:population
* **Genomic Data Source:** GVF (version 1.07)
* **See Also:** <a href="ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/README.populations">1000 Genomes project population codes</a>
* **In Model Version:** 1
* **Example:**
<pre>
    {
        "@context" : "https://www.codamono.com/jsonld/gvf-x1.json",
        "population" : "YRI"
    }
</pre>

### Key: python-callback
* **Description:** Package and module name used for the Python API.
* **Appears in:** context object
* **JSON Type:** string (&ldquo;null&rdquo;, if the Python API was not used)
* **JSON-LD Type:** gfvo-squared:pythonCallback
* **In Model Version:** 1
* **Example:**
<pre>
    {
        "@context" : "https://www.codamono.com/jsonld/vcf-c1.json",
        "python-callback" : "simplepy.simple"
    }
</pre>

### Key: read-pair-span
* **Description:** Global assignment that applies to all features in the data set.
* **Appears in:** meta object
* **JSON Type:** array of numbers
* **JSON-LD Type:** gfvo-squared:readPairSpan
* **Genomic Data Source:** GVF (version 1.07)
* **In Model Version:** 1
* **Example:**
<pre>
    {
        "@context" : "https://www.codamono.com/jsonld/gvf-x1.json",
        "technology-platform" : {
            "SP1" : {
                "read-pair-span" : [
                    135,
                    440
                ]
            }
        }
    }
</pre>

### Key: reads-with-zero-mapping-quality
* **Description:** Number of reads with mapping quality being equal to zero.
* **Appears in:** feature object
* **JSON Type:** number
* **JSON-LD Type:** gfvo-squared:readsWithZeroMappingQuality
* **Genomic Data Source:** VCF (version 4.2)
* **In Model Version:** 1
* **Example:**
<pre>
    {
        "@context" : "https://www.codamono.com/jsonld/vcf-f1.json",
        "reads-with-zero-mapping-quality" : 2
    }
</pre>

### Key: reference
* **Description:** Reference sequence information.
* **Appears in:** feature object
* **JSON Type:** object
* **JSON-LD Type:** gfvo-squared:reference
* **Genomic Data Source:** GVF (version 1.07), VCF (version 4.2)
* **In Model Version:** 1
* **Example:**
<pre>
    {
        "@context" : "https://www.codamono.com/jsonld/vcf-f1.json",
        "reference" : {
            "sequence" : "A"
        }
    }
</pre>

### Key: reference-fasta
* **Description:** A link to a FASTA file that contains reference sequences.
* **Appears in:** meta object
* **JSON Type:** string
* **JSON-LD Type:** gfvo-squared:referenceFASTA
* **Genomic Data Source:** GVF (version 1.07), VCF (version 4.2)
* **In Model Version:** 1
* **Example:**
<pre>
    {
        "@context" : "https://www.codamono.com/jsonld/vcf-x1.json",
        "reference-fasta" : "ftp://ftp-mouse.sanger.ac.uk/ref/GRCm38_68.fa"
    }
</pre>

### Key: runtime
* **Description:** Information about runtime statistics of the BioInterchange software.
* **Appears in:** summary object
* **JSON Type:** object
* **JSON-LD Type:** gfvo-squared:runtime
* **In Model Version:** 1
* **Example:**
<pre>
    {
        "@context" : "https://www.codamono.com/jsonld/biointerchange-s1.json",
         "runtime" : {
            "invocation" : "Fri Jul 17 17:00:54 2015",
            "finish" : "Fri Jul 17 17:00:56 2015",
            "lapsed-seconds" : "2"
        }
    }
</pre>

### Key: samples
* **Description:** Sample information
* **Appears in:** feature object
* **JSON Type:** orray of bject
* **JSON-LD Type:** gfvo-squared:samples
* **Genomic Data Source:** VCF (version 4.2)
* **In Model Version:** 1
* **Example:**
<pre>
    {
        "@context" : "https://www.codamono.com/jsonld/vcf-f1.json",
        "samples" : [
            {
                "id" : "129P2",
                "genotype" : {
                    "sequences" : [
                        "G",
                        "A"
                    ],
                    "phased" : false,
                    "alleles" : "AB"
                }
            }
        ]
    }
</pre>

### Key: samples-with-data
* **Description:** Number of samples with data.
* **Appears in:** feature object
* **JSON Type:** number
* **JSON-LD Type:** gfvo-squared:samplesWithData
* **Genomic Data Source:** VCF (version 4.2)
* **In Model Version:** 1
* **Example:**
<pre>
    {
        "@context" : "https://www.codamono.com/jsonld/vcf-f1.json",
        "samples-with-data" : 1
    }
</pre>

### Key: score-method
* **Description:** Information about scoring methods that are used in a data set.
* **Appears in:** meta object
* **JSON Type:** object
* **JSON-LD Type:** gfvo-squared:scoreMethod
* **Genomic Data Source:** GVF (version 1.07)
* **In Model Version:** 1
* **Example:**
<pre>
    {
        "@context" : "https://www.codamono.com/jsonld/gvf-x1.json",
        "score-method" : {
            "global" : {
                "comment" : "Scores are Phred scaled probabilities of an incorrect sequence_alteration call"
            }
        }
    }
</pre>

### Key: sequence
* **Description:** Genomic or proteomic sequence.
* **Appears in:** feature object
* **JSON Type:** string
* **JSON-LD Type:** gfvo-squared:sequence
* **Genomic Data Source:** GFF3 (version 1.21), GVF (version 1.07), VCF (version 4.2)
* **In Model Version:** 1
* **Example:**
<pre>
    {
        "@context" : "https://www.codamono.com/jsonld/gff3-f1.json",
        "id" : "CDS1",
        "locus" : {
            "landmark" : "chr4",
            "start" : 2,
            "end" : 36,
            "strand" : "-"
        },
        "type" : "CDS",
        "sequence" : "gttcattgctgcctgcatgttcattgtctacctcg"
    }
</pre>

### Key: sequences
* **Description:** A list of genomic or proteomic sequences.
* **Appears in:** feature object
* **JSON Type:** object
* **JSON-LD Type:** gfvo-squared:sequences
* **Genomic Data Source:** VCF (version 4.2)
* **In Model Version:** 1
* **Example:**
<pre>
    {
        "@context" : "https://www.codamono.com/jsonld/vcf-f1.json",
        "samples" : [
            {
                "id" : "129P2",
                "genotype" : {
                    "sequences" : [
                        "G",
                        "A"
                    ],
                    "alleles" : "AB"
                }
            }
        ]
    }
</pre>

### Key: sex
* **Description:** Sex of a sequenced individual.
* **Appears in:** meta object
* **JSON Type:** string
* **JSON-LD Type:** gfvo-squared:sex
* **Genomic Data Source:** GVF (version 1.07)
* **In Model Version:** 1
* **Example:**
<pre>
    {
        "@context" : "https://www.codamono.com/jsonld/gvf-x1.json",
        "sex" : "male"
    }
</pre>

### Key: somatic-mutation
* **Description:** Indicates whether a feature is a somatic mutation.
* **Appears in:** feature object
* **JSON Type:** boolean
* **JSON-LD Type:** gfvo-squared:somaticMutation
* **Genomic Data Source:** VCF (version 4.2)
* **In Model Version:** 1
* **Example:**
<pre>
    {
        "@context" : "https://www.codamono.com/jsonld/vcf-f1.json",
        "somatic-mutation" : true
    }
</pre>

### Key: source
* **Description:** Database or algorithm name that is the source of a feature or data set.
* **Appears in:** feature object
* **JSON Type:** string
* **JSON-LD Type:** gfvo-squared:source
* **Genomic Data Source:** GFF3 (version 1.21), GVF (version 1.07)
* **In Model Version:** 1
* **Example:**
<pre>
    {
        "@context" : "https://www.codamono.com/jsonld/gvf-f1.json",
        "source" : "samtools"
    }
</pre>

### Key: source-method
* **Description:** Information about data sources in a data set.
* **Appears in:** meta object
* **JSON Type:** object
* **JSON-LD Type:** gfvo-squared:sourceMethod
* **Genomic Data Source:** GVF (version 1.07)
* **In Model Version:** 1
* **Example:**
<pre>
    {
        "@context" : "https://www.codamono.com/jsonld/gvf-x1.json",
        "source-method" : {
            "SP6" : {
                "comment" : "Short Elongated Alignment Program (SOAP)",
                "sources" : [
                    "SOAP"
                ],
                "types" : [
                    "SNV"
                ],
                "dbxref" : [
                    "PMID:18227114",
                    "PMID:18987735"
                ]
            }
        }
    }
</pre>

### Key: sources
* **Description:** List of data sources.
* **Appears in:** meta object
* **JSON Type:** array of strings
* **JSON-LD Type:** gfvo-squared:sources
* **Genomic Data Source:** GVF (version 1.07)
* **In Model Version:** 1
* **Example:**
<pre>
    {
        "@context" : "https://www.codamono.com/jsonld/gvf-f1.json",
        "technology-platform" : {
            "SP1" : {
                "sources" : [
                    "Genbank", "dbSNP"
                ]
            }
        }
    }
</pre>

### Key: start
* **Description:** Start coordinate of a feature or alignment.
* **Appears in:** feature object
* **JSON Type:** number
* **JSON-LD Type:** gfvo-squared:start
* **Genomic Data Source:** GFF3 (version 1.21), GVF (version 1.07), VCF (version 4.2)
* **In Model Version:** 1
* **Example:**
<pre>
    {
        "@context" : "https://www.codamono.com/jsonld/gff3-f1.json",
        "locus" : {
            "landmark" : "Chr1",
            "start" : 290207,
            "end" : 291002,
            "strand" : "+"
        },
        "alignment" : {
            "id" : "EST23",
            "start" : 1,
            "end" : 21,
            "strand" : null,
            "cigar-string" : "8M3D6M1I6M"
        }
    }
</pre>

### Key: statistics
* **Description:** Summary statistics about the genomic data that was processed.
* **Appears in:** summary object
* **JSON Type:** object
* **JSON-LD Type:** gfvo-squared:statistics
* **In Model Version:** 1
* **Example:**
<pre>
    {
        "@context" : "https://www.codamono.com/jsonld/biointerchange-s1.json",
        "statistics" : {
            "meta-lines" : 19,
            "meta-lines-filtered" : false,
            "features" : 12,
            "features-filtered" : 0,
            "comment-lines" : 0
        }
    }
</pre>

### Key: strand
* **Description:** Strand on which a feature is located.
* **Appears in:** feature object
* **JSON Type:** string ("+": forward strand, "-": reverse strand)
* **JSON-LD Type:** gfvo-squared:strand
* **Genomic Data Source:** GFF3 (version 1.21), GVF (version 1.07), VCF (version 4.2)
* **In Model Version:** 1
* **Example:**
<pre>
    {
        "@context" : "https://www.codamono.com/jsonld/gff3-f1.json",
        "locus" : {
            "landmark" : "Chr1",
            "start" : 290207,
            "end" : 291002,
            "strand" : "+"
        },
        "alignment" : {
            "id" : "EST23",
            "start" : 1,
            "end" : 21,
            "strand" : "+",
            "cigar-string" : "8M3D6M1I6M"
        }
    }
</pre>

### Key: strand-bias
* **Description:** Strand bias.
* **Appears in:** feature object
* **JSON Type:** number
* **JSON-LD Type:** gfvo-squared:strandBias
* **Genomic Data Source:** VCF (version 4.2)
* **In Model Version:** 1
* **Example:**
<pre>
    {
        "@context" : "https://www.codamono.com/jsonld/vcf-f1.json",
        "strand-bias" : 0.5
    }
</pre>

### Key: technology-platform
* **Description:** Information about the technology platform that was used to create a data set.
* **Appears in:** meta object
* **JSON Type:** object
* **JSON-LD Type:** gfvo-squared:technologyPlatform
* **Genomic Data Source:** GVF (version 1.07)
* **In Model Version:** 1
* **Example:**
<pre>
    {
        "@context" : "https://www.codamono.com/jsonld/gvf-x1.json",
        "technology-platform" : {
            "SP1" : {
                "types" : [
                    "SNV"
                ],
                "comment" : "SNV data described in Wiki."
            }
        }
    }
</pre>

### Key: type
* **Description:** Feature type.
* **Appears in:** feature object
* **JSON Type:** string
* **JSON-LD Type:** gfvo-squared:type
* **Genomic Data Source:** GFF3 (version 1.21), GVF (version 1.07)
* **In Model Version:** 1
* **Example:**
<pre>
    {
        "@context" : "https://www.codamono.com/jsonld/gvf-f1.json",
        "type" : "SNV"
    }
</pre>

### Key: types
* **Description:** List of feature types.
* **Appears in:** meta object
* **JSON Type:** array of strings
* **JSON-LD Type:** gfvo-squared:types
* **Genomic Data Source:** GVF (version 1.07)
* **In Model Version:** 1
* **Example:**
<pre>
    {
        "@context" : "https://www.codamono.com/jsonld/gvf-f1.json",
        "technology-platform" : {
            "SP1" : {
                "types" : [
                    "SNV"
                ],
                "comment" : "SNV data described in Wiki."
            }
        }
    }
</pre>

### Key: user-defined
* **Description:** Key/value pairs that are not defined in the GFF3-, GVF-, or VCF-specification.
* **Appears in:** meta object, feature object
* **JSON Type:** object containing key/value pairs
* **JSON-LD Type:** gfvo-squared:user-defined
* **Genomic Data Source:** GFF3 (version 1.21), GVF (version 1.07), VCF (version 4.2)
* **See Also:** Appendix: unknown keys
* **In Model Version:** 1
* **Example:**
<pre>
    {
        "@context" : "https://www.codamono.com/jsonld/vcf-f1.json",
        "user-defined" : {
            "VDB" : 0.0006,
            "INDEL" : true,
            "DP4" : "0,0,71,0"
        }
    }
</pre>

### Key: user-parameter
* **Description:** &ldquo;-u&rdquo; parameter value when running the BioInterchange software.
* **Appears in:** context object
* **JSON Type:** string (&ldquo;null&rdquo; when no &ldquo;-u&rdquo; was used)
* **JSON-LD Type:** gfvo-squared:user-parameter
* **See Also:** Appendix: unknown keys
* **In Model Version:** 1
* **Example:**
<pre>
    {
        "@context" : "https://www.codamono.com/jsonld/gff3-c1.json",
        "user-parameter" : "Preliminary data."
    }
</pre>

### Key: variants
* **Description:** A collection of sequence variants.
* **Appears in:** feature object
* **JSON Type:** object
* **JSON-LD Type:** gfvo-squared:variants
* **Genomic Data Source:** GVF (version 1.07), VCF (version 4.2)
* **In Model Version:** 1
* **Example:**
<pre>
    {
        "@context" : "https://www.codamono.com/jsonld/gvf-f1.json",
        "variants" : {
            "B" : {
                "sequence" : "T"
            }
        }
    }
</pre>

### Key: vcf-version
* **Description:** Genomics file-format versioning of VCF data sources.
* **Appears in:** meta object
* **JSON Type:** string
* **JSON-LD Type:** gfvo-squared:vcfVersion
* **Genomic Data Source:** VCF (version 4.2)
* **In Model Version:** 1
* **Example:**
<pre>
    {
        "@context" : "https://www.codamono.com/jsonld/vcf-x1.json",
        "vcf-version" : "4.2"
    }
</pre>

### Appendix: unknown keys
* **Description:** Keys that are user defined.
* **Appears in:** feature object, meta object
* **JSON Type:** string
* **JSON-LD Type:** gfvo-squared:unknownProperty-\*
* **Genomic Data Source:** GFF3 (version 1.21), GVF (version 1.07), VCF (version 4.2)
* **In Model Version:** 1
* **Note:** The GFVO<sup>2</sup> type expands to include the unknown key. For example, "VDB" becomes "gfvo-squared:unknownProperty-VDB".
* **Example:**
<pre>
    {
        "@context" : "https://www.codamono.com/jsonld/vcf-f1.json",
        "user-defined" : {
            "VDB" : 0.0006,
            "INDEL" : true,
            "DP4" : "0,0,71,0"
        }
    }
</pre>


