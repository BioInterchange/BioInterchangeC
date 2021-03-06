## Table of Contents

1. [Usage](#usage)
   1. [GFF3, GVF, VCF to JSON](#gff3-gvf-vcf-to-json)
   2. [JSON to GFF3, GVF, VCF](#json-to-gff3-gvf-vcf)
   3. [Python API](#python-api)
   4. [RethinkDB](#rethinkdb)
   5. [MongoDB](#mongodb)
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

BioInterchange is open source now and free to use. If your version of BioInterchange complains about a missing licensing file, then update your BioInterchange installation to the latest verion.

<div class="row">
    <div class="col-xs-12 col-sm-10 col-md-10 col-lg-10 col-xs-offset-0 col-sm-offset-1 col-md-offset-1 col-lg-offset-1">
        <div class="panel panel-info" role="alert">
            <div class="panel-heading">
                Note
            </div>
            <div class="panel-body">
                Trial licenses and licenses for reviewers in the peer-review process are free.
            </div>
        </div>
    </div>
</div>

##### Abbreviations

<dl class="dl-horizontal">
    <dt>GFF3</dt><dd><a href="http://www.sequenceontology.org/gff3.shtml">Generic Feature Format Version 3</a></dd>
    <dt>GVF</dt><dd><a href="http://www.sequenceontology.org/resources/gvf.html">Genome Variation Format</a></dd>
    <dt>JSON</dt><dd><a href="http://json.org/">JavaScript Object Notation</a></dd>
    <dt>JSON-LD</dt><dd><a href="http://json-ld.org/">JSON Linked Data</a></dd>
    <dt>LDJ/LDJSON</dt><dd><a href="https://en.wikipedia.org/wiki/Line_Delimited_JSON">Line Delimited JSON</a></dd>
    <dt>VCF</dt><dd><a href="https://samtools.github.io/hts-specs/VCFv4.2.pdf">Variant Call Format</a></dd>
</dl>

### GFF3, GVF, VCF to JSON

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

### JSON to GFF3, GVF, VCF

Converting JSON objects back to their original genomic file format:

    biointerchange -o example.vcf example.ldj

<div class="row">
    <div class="col-xs-12 col-sm-10 col-md-10 col-lg-10 col-xs-offset-0 col-sm-offset-1 col-md-offset-1 col-lg-offset-1">
        <div class="panel panel-info" role="alert">
            <div class="panel-heading">
                Note
            </div>
            <div class="panel-body">
                Modified JSON objects (for example, when using the Python API) will only translate back correctly when non-standard key/value pairs are put under the &ldquo;user-defined&rdquo; key.
            </div>
        </div>
    </div>
</div>

### Python API

Genomics data can be accessed and modified via the Python API. Each JSON object will be passed on to a Python function and changes made by the Python code will be preserved:

    biointerchange -p simplepy.simple test-data/playground.gff3

<div class="row">
    <div class="col-xs-12 col-sm-10 col-md-10 col-lg-10 col-xs-offset-0 col-sm-offset-1 col-md-offset-1 col-lg-offset-1">
        <div class="panel panel-info" role="alert">
            <div class="panel-heading">
                Compatibility
            </div>
            <div class="panel-body">
                Requires Python 3.4 or newer.
            </div>
        </div>
    </div>
</div>

<dl>
    <dt>Example: Directory structure.</dt>
    <dd markdown="1">
* ./simplepy/\_\_init\_\_.py
* ./simplepy/simple.py
</dd>
    <dt>Example: Environment variables.</dt>
    <dd markdown="1">
* PYTHONHOME needs to be set to where Python is installed on your system
* PYTHONPATH needs to include "." for this example, so that "simple.py" can be found
</dd>
    <dt>Example: &ldquo;simple.py&rdquo; source code.</dt>
    <dd markdown="1">
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
</dd>
</dl>

### RethinkDB

Take the cat lovers example from the main page:

    wget ftp://ftp.ensembl.org/pub/release-81/variation/vcf/felis_catus/Felis_catus_incl_consequences.vcf.gz
    gunzip Felis_catus_incl_consequences.vcf.gz
    biointerchange -o Felis_catus_incl_consequences.ldj Felis_catus_incl_consequences.vcf

Import line-delimited JSON-LD documents into RethinkDB:

    rethinkdb import -f Felis_catus_incl_consequences.ldj --table genomics.felis_catus

Check that the data is actually in the database using RethinkDB's "Data Explorer" (for local installation at `http://localhost:8080/#dataexplorer`):

    r.db('genomics').table('felis_catus')

Voila!

### MongoDB

<div class="row">
    <div class="col-xs-12 col-sm-10 col-md-10 col-lg-10 col-xs-offset-0 col-sm-offset-1 col-md-offset-1 col-lg-offset-1">
        <div class="panel panel-info" role="alert">
            <div class="panel-heading">
                Note
            </div>
            <div class="panel-body">
                MongoDB can be a bit much to begin with. Have a look at RethinkDB for starters! RethinkDB has a great intuitive web-interface, easy-peasy indexing and sharding support, and a cool query language (ReQL).
            </div>
        </div>
    </div>
</div>

Take the cat lovers example from the main page:

    wget ftp://ftp.ensembl.org/pub/release-81/variation/vcf/felis_catus/Felis_catus_incl_consequences.vcf.gz
    gunzip Felis_catus_incl_consequences.vcf.gz
    biointerchange -o Felis_catus_incl_consequences.ldj Felis_catus_incl_consequences.vcf

Import line-delimited JSON-LD documents into MongoDB:

    mongoimport --db genomics --collection felis_catus --type json --file Felis_catus_incl_consequences.ldj

Check that the data is actually in the database:

    mongo genomics
    > db.felis_catus.find()

Voila!

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

<div class="row">
    <div class="col-xs-12 col-sm-10 col-md-10 col-lg-10 col-xs-offset-0 col-sm-offset-1 col-md-offset-1 col-lg-offset-1">
        <div class="panel panel-info" role="alert">
            <div class="panel-heading">
                Note
            </div>
            <div class="panel-body">
                In a nutshell: meta objects and feature objects contain the genomics data from GFF3, GVF, and VCF files; these objects are referencing each other (they are linked). Context and summary objects stand on their own.
            </div>
        </div>
    </div>
</div>

**All** objects contain a "@context" key and a "@type" key, which are called context key and type key in the following.

The context key turns the JSON objects into JSON-LD objects. Never heard of JSON-LD? Just skip to the next part and ignore the "@context" key. This is what makes JSON-LD so great: you can handle JSON-LD objects just like JSON objects.

If you do want to make use of the context key, then feed the JSON-LD objects to a Linked Data tool and it will annotate key/value pairs with type information. JSON-LD objects can also be turned into Triple Store compatible data formats, such as RDF N-Triples, RDF N-Quads, and RDF Turtle. Want to see a real-life example of this magic? Head over to the [JSON-LD Playground](http://json-ld.org/playground/) and copy/paste any JSON-LD object into the "JSON-LD Input" text field: the "N-Quads" tab will instantly show a triples representation of the JSON-LD object that can be loaded into a triple store.

<div class="row">
    <div class="col-xs-12 col-sm-10 col-md-10 col-lg-10 col-xs-offset-0 col-sm-offset-1 col-md-offset-1 col-lg-offset-1">
        <div class="panel panel-info" role="alert">
            <div class="panel-heading">
                Note
            </div>
            <div class="panel-body">
                All JSON-LD annotations are making use of GFVO<sup>2</sup>. Is it legit? Well, it is based on all the lessons learned when designing the <a href="https://peerj.com/articles/933/">Genomic Feature and Variation Ontology</a> (<a href="https://peerj.com/articles/933/">GFVO</a>). So, there is a team behind GFVO<sup>2</sup> too? Yes, stay tuned&hellip;
            </div>
        </div>
    </div>
</div>

#### Context Objects

Context objects tell you something about the environment in which the following JSON objects ("meta objects", "feature objects" and "summary objects") were created.

In essence, a context objects contains information about:

* which version of BioInterchange was used
* the input filename
* the filetype of the input (GFF3, GVF, VCF)
* what the output filename was (unless output went to the console)
* whether the Python API was utilized (if so, name of the Python module)
* additional user-defined parameters

<dl>
    <dt>Example</dt>
    <dd markdown="1">
    {
        "@context" : "https://www.codamono.com/jsonld/vcf-c1.json",
        "biointerchange-version" : "2.0.0+36",
        "input-filetype" : "VCF",
        "output-file" : null,
        "input-file" : "test-data/playground-vcf.vcf",
        "python-callback" : null,
        "user-defined" : null
    }
</dd>
</dl>

Detailed information about each key/value-pair can be found in the JSON Reference Card section.

#### Meta Objects

For each genomic file, one meta object is being created. The meta object contain data of information/pragma lines.

Meta objects are helpful for determining data provenance, establishing links to ontologies that were used, and providing extra annotations that were not attached to features to reduce data redundancy.

For example, the following example provides a textual descriptions for analytic filters that have been applied in a VCF file. Instead of adding this description to every genomic feature for which the filter applies, the textual description is only given once in this meta object, where it can be looked up via the filter keys ("MinAB", "MinDP", "MinMQ", and "Qual").

<dl>
    <dt>Example:</dt>
    <dd markdown="1">
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
</dd>
</dl>

#### Feature Objects

Sequences, sequence variations, sequence annotations, genotyping samples, etc., are represented by feature objects.

##### Basic Feature Information

Most basic information about features includes an identifier, a genomic locus, provenane ("SGRP" stands for the Saccharomyces Genome Resequencing Project), feature type information ("SNV" -- single nucleotide variant; a Sequence Ontology term), and references to external databases (SGRP, again, and European Molecular Biology Laboratory accession).

<dl>
    <dt>Example</dt>
    <dd markdown="1">
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
</dd>
</dl>

##### Variant and Reference Sequences

Information about reference sequences is stored under the "reference" key. Variations are stored under the "variants" key and allele specific informations is labeled by "B", "C", etc. ("A" is denoting the reference, which is not explicitly labeled).

<dl>
    <dt>Example: GVF (basic information omitted)</dt>
    <dd markdown="1">
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
</dd>
    <dt>Example: VCF (basic information omitted)</dt>
    <dd markdown="1">
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
</dd>
</dl>

##### Genomic/Genotyping Samples

VCF genomic files contain information about samples. Reference, variant, and other information is repeated for each sample.

<dl>
    <dt>Example: VCF (basic and non-sample specific information omitted)</dt>
    <dd markdown="1">
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
</dd>
</dl>

#### Summary Objects

Summary objects wrap-up everything that came beforehand. Namely:

* statistics about the genomics data
* runtime information

Genomics data statistics capture how many comments were seen, how much metadata was read, and how much actual genomic features were processed.

Runtime information tell you when BioInterchange was invoked, when it finished processing the data, and how many seconds it took to process the data.

<dl>
    <dt>Example:</dt>
    <dd markdown="1">
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
</dd>
</dl>

## Python API Reference Cards

<div class="panel panel-default">
    <div class="panel-heading">Function: setup(context, meta)</div>
    <div class="panel-body">
        <dl class="dl-horizontal">
            <dt>Description</dt>
            <dd>Initialization function.</dd>
            <dt>Parameter &ldquo;context&rdquo;</dt>
            <dd>Context object as a Python Dict-instance. Modifications to this dict will have no effect on the output of BioInterchange.</dd>
            <dt>Parameter &ldquo;meta&rdquo;</dt>
            <dd>Meta object as a Python Dict-instance.</dd>
            <dt>Returns</dt>
            <dd>The original meta object or a modified version of it. If None is returned, then the output of the meta object will be suppressed by BioInterchange</dd>
        </dl>
        <dl>
            <dt>Example:</dt>
            <dd markdown="1">
~~~ python
def setup(context, meta):
    return meta
~~~
</dd>
        </dl>
    </div>
</div>

<div class="panel panel-default">
    <div class="panel-heading">Function: cleanup(summary)</div>
    <div class="panel-body">
        <dl class="dl-horizontal">
            <dt>Description</dt>
            <dd>Cleanup and finalization function.</dd>
            <dt>Parameter &ldquo;summary&rdquo;</dt>
            <dd>Summary object as a Python Dict-instance.</dd>
            <dt>Returns</dt>
            <dd>The original summary object or a modified version of it. If None is returned, then the output of the summary object will be suppressed by BioInterchange</dd>
        </dl>
        <dl>
            <dt>Example:</dt>
            <dd markdown="1">
~~~ python
def cleanup(summary):
    return summary
~~~
</dd>
        </dl>
    </div>
</div>

<div class="panel panel-default">
    <div class="panel-heading">Function: process_feature(feature)</div>
    <div class="panel-body">
        <dl class="dl-horizontal">
            <dt>Description</dt>
            <dd>Genomic data processing function.</dd>
            <dt>Parameter &ldquo;featurey&rdquo;</dt>
            <dd>Feature object as a Python Dict-instance.</dd>
            <dt>Returns</dt>
            <dd>The original feature object or a modified version of it. If None is returned, then the output of the feature object will be suppressed by BioInterchange</dd>
        </dl>
        <dl>
            <dt>Example:</dt>
            <dd markdown="1">
~~~ python
def process_feature(feature):
    return feature
~~~
</dd>
        </dl>
    </div>
</div>

## JSON Reference Cards

<div class="panel panel-default">
    <div class="panel-heading">Key: B, C, etc.</div>
    <div class="panel-body">
        <dl class="dl-horizontal">
            <dt>Description</dt>
            <dd>Alternative allele information (reference is implicitly labeled "A").</dd>
            <dt>Appears in</dt>
            <dd>feature object</dd>
            <dt>JSON Type</dt>
            <dd>object</dd>
            <dt>JSON-LD Type</dt>
            <dd>gfvo-squared:alleleB, gfvo-squared:alleleC, etc.</dd>
            <dt>Genomic Data Source</dt>
            <dd>GVF (version 1.07), VCF (version 4.2)</dd>
            <dt>In Model Version</dt>
            <dd>1</dd>
        </dl>
        <dl>
            <dt>Example: Variant information on a feature level.</dt>
            <dd markdown="1">
    {
        "@context" : "https://www.codamono.com/jsonld/vcf-f1.json",
        "variants" : {
            "B" : {
                "sequence" : "G"
            }
        }
    }
</dd>
            <dt>Example: Variant information on a sample level.</dt>
            <dd markdown="1">
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
</dd>
        </dl>
    </div>
</div>

<div class="panel panel-default">
    <div class="panel-heading">Key: AA, AB, BB, AC, etc.</div>
    <div class="panel-body">
        <dl class="dl-horizontal">
            <dt>Description</dt>
            <dd>Generic genotype information.</dd>
            <dt>Appears in</dt>
            <dd>feature object</dd>
            <dt>JSON Type</dt>
            <dd>object</dd>
            <dt>JSON-LD Type</dt>
            <dd>gfvo-squared:genotypeAA, gfvo-squared:genotypeAB, etc.</dd>
            <dt>Genomic Data Source</dt>
            <dd>GVF (version 1.07), VCF (version 4.2)</dd>
            <dt>See Also</dt>
            <dd>genotype</dd>
            <dt>In Model Version</dt>
            <dd>1</dd>
        </dl>
        <dl>
            <dt>Example</dt>
            <dd markdown="1">
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
</dd>
        </dl>
    </div>
</div>

<div class="panel panel-default">
    <div class="panel-heading">Key: affected-features</div>
    <div class="panel-body">
        <dl class="dl-horizontal">
            <dt>Description</dt>
            <dd>Identifiers of features that are affected by a variant effect.</dd>
            <dt>Appears in</dt>
            <dd>feature object</dd>
            <dt>JSON Type</dt>
            <dd>array of strings</dd>
            <dt>JSON-LD Type</dt>
            <dd>gfvo-squared:affectedFeatures</dd>
            <dt>Genomic Data Source</dt>
            <dd>GVF (version 1.07)</dd>
            <dt>In Model Version</dt>
            <dd>1</dd>
        </dl>
        <dl>
            <dt>Example</dt>
            <dd markdown="1">
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
</dd>
        </dl>
    </div>
</div>

<div class="panel panel-default">
    <div class="panel-heading">Key: affected-feature-type</div>
    <div class="panel-body">
        <dl class="dl-horizontal">
            <dt>Description</dt>
            <dd>Type of the features that are affected by a variant effect.</dd>
            <dt>Appears in</dt>
            <dd>feature object</dd>
            <dt>JSON Type</dt>
            <dd>string</dd>
            <dt>JSON-LD Type</dt>
            <dd>gfvo-squared:affectedFeatureType</dd>
            <dt>Genomic Data Source</dt>
            <dd>GVF (version 1.07)</dd>
            <dt>In Model Version</dt>
            <dd>1</dd>
        </dl>
        <dl>
            <dt>Example</dt>
            <dd markdown="1">
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
</dd>
        </dl>
    </div>
</div>

<div class="panel panel-default">
    <div class="panel-heading">Key: alias</div>
    <div class="panel-body">
        <dl class="dl-horizontal">
            <dt>Description</dt>
            <dd>Aliases of a feature.</dd>
            <dt>Appears in</dt>
            <dd>feature object</dd>
            <dt>JSON Type</dt>
            <dd>array of strings</dd>
            <dt>JSON-LD Type</dt>
            <dd>gfvo-squared:alias</dd>
            <dt>Genomic Data Source</dt>
            <dd>GFF3 (version 1.21), GVF (version 1.07), VCF (version 4.2)</dd>
            <dt>Note</dt>
            <dd>In VCF, only the first entry of the "ID" column becomes an "id" in JSON-LD, whereas second, third, etc., entries are interpreted as "alias" in JSON-LD.</dd>
            <dt>In Model Version</dt>
            <dd>1</dd>
        </dl>
        <dl>
            <dt>Example</dt>
            <dd markdown="1">
    {
        "@context" : "https://www.codamono.com/jsonld/gff3-f1.json",
        "id" : "rs123",
        "alias" : [
            "feat12",
            "feat12-1"
        ]
    }
</dd>
        </dl>
    </div>
</div>

<div class="panel panel-default">
    <div class="panel-heading">Key: alignment</div>
    <div class="panel-body">
        <dl class="dl-horizontal">
            <dt>Description</dt>
            <dd>A sequence alignment.</dd>
            <dt>Appears in</dt>
            <dd>feature object</dd>
            <dt>JSON Type</dt>
            <dd>object</dd>
            <dt>JSON-LD Type</dt>
            <dd>gfvo-squared:alignment</dd>
            <dt>Genomic Data Source</dt>
            <dd>GFF3 (version 1.21), GVF (version 1.07), VCF (version 4.2)</dd>
            <dt>In Model Version</dt>
            <dd>1</dd>
        </dl>
        <dl>
            <dt>Example</dt>
            <dd markdown="1">
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
</dd>
        </dl>
    </div>
</div>

<div class="panel panel-default">
    <div class="panel-heading">Key: allele-count</div>
    <div class="panel-body">
        <dl class="dl-horizontal">
            <dt>Description</dt>
            <dd>Number of alleles in a genotype.</dd>
            <dt>Appears in</dt>
            <dd>feature object</dd>
            <dt>JSON Type</dt>
            <dd>number</dd>
            <dt>JSON-LD Type</dt>
            <dd>gfvo-squared:alleleCount</dd>
            <dt>Genomic Data Source</dt>
            <dd>VCF (version 4.2)</dd>
            <dt>In Model Version</dt>
            <dd>1</dd>
        </dl>
        <dl>
            <dt>Example</dt>
            <dd markdown="1">
    {
        "@context" : "https://www.codamono.com/jsonld/vcf-f1.json",
        "variants" : {
            "B" : {
                "sequence" : "GA",
                "allele-count" : 36
            }
        }
    }
</dd>
        </dl>
    </div>
</div>

<div class="panel panel-default">
    <div class="panel-heading">Key: allele-count-expected</div>
    <div class="panel-body">
        <dl class="dl-horizontal">
            <dt>Description</dt>
            <dd>Expected alternate allele count.</dd>
            <dt>Appears in</dt>
            <dd>feature object</dd>
            <dt>JSON Type</dt>
            <dd>number</dd>
            <dt>JSON-LD Type</dt>
            <dd>gfvo-squared:alleleCountExpected</dd>
            <dt>Genomic Data Source</dt>
            <dd>VCF (version 4.2)</dd>
            <dt>In Model Version</dt>
            <dd>1</dd>
        </dl>
        <dl>
            <dt>Example</dt>
            <dd markdown="1">
    {
        "@context" : "https://www.codamono.com/jsonld/vcf-f1.json",
        "variants" : {
            "B" : {
                "sequence" : "GA",
                "allele-count-expected" : 4
            }
        }
    }
</dd>
        </dl>
    </div>
</div>

<div class="panel panel-default">
    <div class="panel-heading">Key: allele-format</div>
    <div class="panel-body">
        <dl class="dl-horizontal">
            <dt>Description</dt>
            <dd>Information about key/values used to describe alleles in a data set.</dd>
            <dt>Appears in</dt>
            <dd>meta object</dd>
            <dt>JSON Type</dt>
            <dd>object</dd>
            <dt>JSON-LD Type</dt>
            <dd>gfvo-squared:alleleFormat</dd>
            <dt>Genomic Data Source</dt>
            <dd>VCF (version 4.2)</dd>
            <dt>In Model Version</dt>
            <dd>1</dd>
        </dl>
        <dl>
            <dt>Example</dt>
            <dd markdown="1">
    {
        "@context" : "https://www.codamono.com/jsonld/vcf-x1.json",
        "allele-format" : {
            "DEL" : {
                "description" : "Deletion"
            }
        }
    }
</dd>
        </dl>
    </div>
</div>

<div class="panel panel-default">
    <div class="panel-heading">Key: allele-frequency</div>
    <div class="panel-body">
        <dl class="dl-horizontal">
            <dt>Description</dt>
            <dd>Frequency of an allele in a genotype.</dd>
            <dt>Appears in</dt>
            <dd>feature object</dd>
            <dt>JSON Type</dt>
            <dd>number</dd>
            <dt>JSON-LD Type</dt>
            <dd>gfvo-squared:alleleFrequency</dd>
            <dt>Genomic Data Source</dt>
            <dd>VCF (version 4.2)</dd>
            <dt>In Model Version</dt>
            <dd>1</dd>
        </dl>
        <dl>
            <dt>Example</dt>
            <dd markdown="1">
    {
        "@context" : "https://www.codamono.com/jsonld/vcf-f1.json",
        "variants" : {
            "B" : {
                "sequence" : "GA",
                "allele-frequency" : 1
            }
        }
    }
</dd>
        </dl>
    </div>
</div>

<div class="panel panel-default">
    <div class="panel-heading">Key: allele-total-number</div>
    <div class="panel-body">
        <dl class="dl-horizontal">
            <dt>Description</dt>
            <dd>Total number of alleles.</dd>
            <dt>Appears in</dt>
            <dd>feature object</dd>
            <dt>JSON Type</dt>
            <dd>number</dd>
            <dt>JSON-LD Type</dt>
            <dd>gfvo-squared:alleleTotalNumber</dd>
            <dt>Genomic Data Source</dt>
            <dd>VCF (version 4.2)</dd>
            <dt>In Model Version</dt>
            <dd>1</dd>
        </dl>
        <dl>
            <dt>Example</dt>
            <dd markdown="1">
    {
        "@context" : "https://www.codamono.com/jsonld/vcf-f1.json",
        "allele-total-number" : 36
    }
</dd>
        </dl>
    </div>
</div>

<div class="panel panel-default">
    <div class="panel-heading">Key: ancestral-allele</div>
    <div class="panel-body">
        <dl class="dl-horizontal">
            <dt>Description</dt>
            <dd>Sequence of an ancestral allele.</dd>
            <dt>Appears in</dt>
            <dd>feature object</dd>
            <dt>JSON Type</dt>
            <dd>string</dd>
            <dt>JSON-LD Type</dt>
            <dd>gfvo-squared:ancestralAllele</dd>
            <dt>Genomic Data Source</dt>
            <dd>VCF (version 4.2)</dd>
            <dt>In Model Version</dt>
            <dd>1</dd>
        </dl>
        <dl>
            <dt>Example</dt>
            <dd markdown="1">
    {
        "@context" : "https://www.codamono.com/jsonld/vcf-f1.json",
        "ancestral-allele" : "C"
    }
</dd>
        </dl>
    </div>
</div>

<div class="panel panel-default">
    <div class="panel-heading">Key: annotation-format</div>
    <div class="panel-body">
        <dl class="dl-horizontal">
            <dt>Description</dt>
            <dd>Information about keys used to annotate (filter) features in a data set.</dd>
            <dt>Appears in</dt>
            <dd>meta object</dd>
            <dt>JSON Type</dt>
            <dd>object</dd>
            <dt>JSON-LD Type</dt>
            <dd>gfvo-squared:annotationFormat</dd>
            <dt>Genomic Data Source</dt>
            <dd>VCF (version 4.2)</dd>
            <dt>In Model Version</dt>
            <dd>1</dd>
        </dl>
        <dl>
            <dt>Example</dt>
            <dd markdown="1">
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
</dd>
        </dl>
    </div>
</div>

<div class="panel panel-default">
    <div class="panel-heading">Key: annotations</div>
    <div class="panel-body">
        <dl class="dl-horizontal">
            <dt>Description</dt>
            <dd>Tags (or labels) assigned to a feature.</dd>
            <dt>Appears in</dt>
            <dd>feature object</dd>
            <dt>JSON Type</dt>
            <dd>array of strings</dd>
            <dt>JSON-LD Type</dt>
            <dd>gfvo-squared:annotations</dd>
            <dt>Genomic Data Source</dt>
            <dd>GVF (version 1.07), VCF (version 4.2)</dd>
            <dt>Note</dt>
            <dd>Annotations are represented by structured pragma statements in GVF and "FILTER" information in VCF.</dd>
            <dt>In Model Version</dt>
            <dd>1</dd>
        </dl>
        <dl>
            <dt>Example: GVF</dt>
            <dd markdown="1">
    {
        "@context" : "https://www.codamono.com/jsonld/gvf-f1.json",
        "annotations" : [
            "SP1"
        ]
    }
</dd>
            <dt>Annotation Source Example (Meta Document)</dt>
            <dd markdown="1">
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
</dd>
        </dl>
    </div>
</div>

<div class="panel panel-default">
    <div class="panel-heading">Key: attribute-method</div>
    <div class="panel-body">
        <dl class="dl-horizontal">
            <dt>Description</dt>
            <dd>Information about attributes in the data set.</dd>
            <dt>Appears in</dt>
            <dd>meta object</dd>
            <dt>JSON Type</dt>
            <dd>object</dd>
            <dt>JSON-LD Type</dt>
            <dd>gfvo-squared:attributeMethod</dd>
            <dt>Genomic Data Source</dt>
            <dd>GVF (version 1.07)</dd>
            <dt>In Model Version</dt>
            <dd>1</dd>
        </dl>
        <dl>
            <dt>Example</dt>
            <dd markdown="1">
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
</dd>
        </dl>
    </div>
</div>

<div class="panel panel-default">
    <div class="panel-heading">Key: average-coverage</div>
    <div class="panel-body">
        <dl class="dl-horizontal">
            <dt>Description</dt>
            <dd>Average read coverage.</dd>
            <dt>Appears in</dt>
            <dd>meta object</dd>
            <dt>JSON Type</dt>
            <dd>number</dd>
            <dt>JSON-LD Type</dt>
            <dd>gfvo-squared:averageCoverage</dd>
            <dt>Genomic Data Source</dt>
            <dd>GVF (version 1.07)</dd>
            <dt>In Model Version</dt>
            <dd>1</dd>
        </dl>
        <dl>
            <dt>Example</dt>
            <dd markdown="1">
    {
        "@context" : "https://www.codamono.com/jsonld/gvf-x1.json",
        "technology-platform" : {
            "SP1" : {
                "average-coverage" : 36
            }
        }
    }
</dd>
        </dl>
    </div>
</div>

<div class="panel panel-default">
    <div class="panel-heading">Key: base-quality-rms</div>
    <div class="panel-body">
        <dl class="dl-horizontal">
            <dt>Description</dt>
            <dd>Root-mean-square base quality of a genomic position.</dd>
            <dt>Appears in</dt>
            <dd>feature object</dd>
            <dt>JSON Type</dt>
            <dd>number</dd>
            <dt>JSON-LD Type</dt>
            <dd>gfvo-squared:baseQualityRMS</dd>
            <dt>Genomic Data Source</dt>
            <dd>VCF (version 4.2)</dd>
            <dt>In Model Version</dt>
            <dd>1</dd>
        </dl>
        <dl>
            <dt>Example</dt>
            <dd markdown="1">
    {
        "@context" : "https://www.codamono.com/jsonld/vcf-f1.json",
        "base-quality-rms" : 3
    }
</dd>
        </dl>
    </div>
</div>

<div class="panel panel-default">
    <div class="panel-heading">Key: biointerchange-version</div>
    <div class="panel-body">
        <dl class="dl-horizontal">
            <dt>Description</dt>
            <dd><a href="http://semver.org/">Semantic Version</a> number of the BioInterchange software that was used to create the data set.</dd>
            <dt>Appears in</dt>
            <dd>context object</dd>
            <dt>JSON Type</dt>
            <dd>string</dd>
            <dt>JSON-LD Type</dt>
            <dd>gfvo-squared:biointerchangeVersion</dd>
            <dt>In Model Version</dt>
            <dd>1</dd>
        </dl>
        <dl>
            <dt>Example</dt>
            <dd markdown="1">
    {
        "@context" : "https://www.codamono.com/jsonld/vcf-c1.json",
        "biointerchange-version" : "2.0.0+36"
    }
</dd>
        </dl>
    </div>
</div>

<div class="panel panel-default">
    <div class="panel-heading">Key: breakpoint-fasta</div>
    <div class="panel-body">
        <dl class="dl-horizontal">
            <dt>Description</dt>
            <dd>A link to a FASTA file that contains sequences specific to breakpoints.</dd>
            <dt>Appears in</dt>
            <dd>meta object</dd>
            <dt>JSON Type</dt>
            <dd>string</dd>
            <dt>JSON-LD Type</dt>
            <dd>gfvo-squared:breakpointFASTA</dd>
            <dt>Genomic Data Source</dt>
            <dd>VCF (version 4.2)</dd>
            <dt>In Model Version</dt>
            <dd>1</dd>
        </dl>
        <dl>
            <dt>Example</dt>
            <dd markdown="1">
    {
        "@context" : "https://www.codamono.com/jsonld/vcf-x1.json",
        "breakpoint-fasta" : "ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/sv/breakpoint_assemblies.fasta"
    }
</dd>
        </dl>
    </div>
</div>

<div class="panel panel-default">
    <div class="panel-heading">Key: build</div>
    <div class="panel-body">
        <dl class="dl-horizontal">
            <dt>Description</dt>
            <dd>Name of a genomic build.</dd>
            <dt>Appears in</dt>
            <dd>meta object</dd>
            <dt>JSON Type</dt>
            <dd>string</dd>
            <dt>JSON-LD Type</dt>
            <dd>gfvo-squared:build</dd>
            <dt>Genomic Data Source</dt>
            <dd>GFF3 (version 1.21), GVF (version 1.07)</dd>
            <dt>In Model Version</dt>
            <dd>1</dd>
        </dl>
        <dl>
            <dt>Example</dt>
            <dd markdown="1">
    {
        "@context" : "https://www.codamono.com/jsonld/gvf-f1.json",
        "genome-build" : [
            {
                "source" : "NCBI",
                "build" : "B36"
            }
        ]
    }
</dd>
        </dl>
    </div>
</div>

<div class="panel panel-default">
    <div class="panel-heading">Key: cigar-string</div>
    <div class="panel-body">
        <dl class="dl-horizontal">
            <dt>Description</dt>
            <dd>A CIGAR formatted sequence alignment string.</dd>
            <dt>Appears in</dt>
            <dd>feature object</dd>
            <dt>JSON Type</dt>
            <dd>string in CIGAR format</dd>
            <dt>JSON-LD Type</dt>
            <dd>gfvo-squared:cigarString</dd>
            <dt>Genomic Data Source</dt>
            <dd>GFF3 (version 1.21), GVF (version 1.07), VCF (version 4.2)</dd>
            <dt>Note</dt>
            <dd>CIGAR strings from GFF3/GVF are reformatted to match the CIGAR standard: integer followed by a character. When translating JSON-LD back to GFF3/GVF, the alternative GFF3/GVF-specific format is substituted again.</dd>
            <dt>In Model Version</dt>
            <dd>1</dd>
        </dl>
        <dl>
            <dt>Example</dt>
            <dd markdown="1">
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
</dd>
        </dl>
    </div>
</div>

<div class="panel panel-default">
    <div class="panel-heading">Key: codon</div>
    <div class="panel-body">
        <dl class="dl-horizontal">
            <dt>Description</dt>
            <dd>A codon sequence.</dd>
            <dt>Appears in</dt>
            <dd>feature object</dd>
            <dt>JSON Type</dt>
            <dd>string</dd>
            <dt>JSON-LD Type</dt>
            <dd>gfvo-squared:codon</dd>
            <dt>Genomic Data Source</dt>
            <dd>GVF (version 1.07)</dd>
            <dt>In Model Version</dt>
            <dd>1</dd>
        </dl>
        <dl>
            <dt>Example</dt>
            <dd markdown="1">
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
</dd>
        </dl>
    </div>
</div>

<div class="panel panel-default">
    <div class="panel-heading">Key: codon-phase</div>
    <div class="panel-body">
        <dl class="dl-horizontal">
            <dt>Description</dt>
            <dd>Phase of a coding sequence.</dd>
            <dt>Appears in</dt>
            <dd>feature object</dd>
            <dt>JSON Type</dt>
            <dd>number (values: 0, 1, 2)</dd>
            <dt>JSON-LD Type</dt>
            <dd>gfvo-squared:codon-phase</dd>
            <dt>Genomic Data Source</dt>
            <dd>GFF3 (version 1.21)</dd>
            <dt>In Model Version</dt>
            <dd>1</dd>
        </dl>
        <dl>
            <dt>Example</dt>
            <dd markdown="1">
    {
        "@context" : "https://www.codamono.com/jsonld/gff3-f1.json",
        "codon-phase" : 0
    }
</dd>
        </dl>
    </div>
</div>

<div class="panel panel-default">
    <div class="panel-heading">Key: comment</div>
    <div class="panel-body">
        <dl class="dl-horizontal">
            <dt>Description</dt>
            <dd>A free-text comment.</dd>
            <dt>Appears in</dt>
            <dd>feature object, meta object</dd>
            <dt>JSON Type</dt>
            <dd>string</dd>
            <dt>JSON-LD Type</dt>
            <dd>gfvo-squared:comment</dd>
            <dt>Genomic Data Source</dt>
            <dd>GFF3 (version 1.21), GVF (version 1.07), VCF (version 4.2)</dd>
            <dt>In Model Version</dt>
            <dd>1</dd>
        </dl>
        <dl>
            <dt>Example</dt>
            <dd markdown="1">
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
</dd>
        </dl>
    </div>
</div>

<div class="panel panel-default">
    <div class="panel-heading">Key: comment-lines</div>
    <div class="panel-body">
        <dl class="dl-horizontal">
            <dt>Description</dt>
            <dd>Number of comment lines that were read from a genomics file.</dd>
            <dt>Appears in</dt>
            <dd>summary object</dd>
            <dt>JSON Type</dt>
            <dd>number</dd>
            <dt>JSON-LD Type</dt>
            <dd>gfvo-squared:commentLines</dd>
            <dt>In Model Version</dt>
            <dd>1</dd>
        </dl>
        <dl>
            <dt>Example:</dt>
            <dd markdown="1">
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
</dd>
        </dl>
    </div>
</div>

<div class="panel panel-default">
    <div class="panel-heading">Key: contig</div>
    <div class="panel-body">
        <dl class="dl-horizontal">
            <dt>Description</dt>
            <dd>Information about a continuous sequence region.</dd>
            <dt>Appears in</dt>
            <dd>meta object</dd>
            <dt>JSON Type</dt>
            <dd>object</dd>
            <dt>JSON-LD Type</dt>
            <dd>gfvo-squared:contig</dd>
            <dt>Genomic Data Source</dt>
            <dd>GFF3 (version 1.21), GVF (version 1.07), VCF (version 4.2)</dd>
            <dt>In Model Version</dt>
            <dd>1</dd>
        </dl>
        <dl>
            <dt>Example</dt>
            <dd markdown="1">
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
</dd>
        </dl>
    </div>
</div>

<div class="panel panel-default">
    <div class="panel-heading">Key: data-source</div>
    <div class="panel-body">
        <dl class="dl-horizontal">
            <dt>Description</dt>
            <dd>Information about the source of genomic/proteomic features.</dd>
            <dt>Appears in</dt>
            <dd>meta object</dd>
            <dt>JSON Type</dt>
            <dd>object</dd>
            <dt>JSON-LD Type</dt>
            <dd>gfvo-squared:dataSource</dd>
            <dt>Genomic Data Source</dt>
            <dd>GVF (version 1.07)</dd>
            <dt>In Model Version</dt>
            <dd>1</dd>
        </dl>
        <dl>
            <dt>Example</dt>
            <dd markdown="1">
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
</dd>
        </dl>
    </div>
</div>

<div class="panel panel-default">
    <div class="panel-heading">Key: data-type</div>
    <div class="panel-body">
        <dl class="dl-horizontal">
            <dt>Description</dt>
            <dd>Type of data that is presented by a data source.</dd>
            <dt>Appears in</dt>
            <dd>meta object</dd>
            <dt>JSON Type</dt>
            <dd>string</dd>
            <dt>JSON-LD Type</dt>
            <dd>gfvo-squared:dataType</dd>
            <dt>Genomic Data Source</dt>
            <dd>GVF (version 1.07)</dd>
            <dt>See Also</dt>
            <dd>Key: data-source</dd>
            <dt>In Model Version</dt>
            <dd>1</dd>
        </dl>
        <dl>
            <dt>Example</dt>
            <dd markdown="1">
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
</dd>
        </dl>
    </div>
</div>

<div class="panel panel-default">
    <div class="panel-heading">Key: dbxref</div>
    <div class="panel-body">
        <dl class="dl-horizontal">
            <dt>Description</dt>
            <dd>External database references.</dd>
            <dt>Appears in</dt>
            <dd>feature object, meta object</dd>
            <dt>JSON Type</dt>
            <dd>string</dd>
            <dt>JSON-LD Type</dt>
            <dd>gfvo-squared:dbxref</dd>
            <dt>Genomic Data Source</dt>
            <dd>GFF3 (version 1.21), GVF (version 1.07)</dd>
            <dt>In Model Version</dt>
            <dd>1</dd>
        </dl>
        <dl>
            <dt>Example</dt>
            <dd markdown="1">
    {
        "@context" : "https://www.codamono.com/jsonld/gvf-f1.json",
        "dbxref" : [
            "SGRP:s01-675",
            "EMBL:AA816246"
        ]
    }
</dd>
        </dl>
    </div>
</div>

<div class="panel panel-default">
    <div class="panel-heading">Key: depth</div>
    <div class="panel-body">
        <dl class="dl-horizontal">
            <dt>Description</dt>
            <dd>Read depth.</dd>
            <dt>Appears in</dt>
            <dd>feature object</dd>
            <dt>JSON Type</dt>
            <dd>number</dd>
            <dt>JSON-LD Type</dt>
            <dd>gfvo-squared:depth</dd>
            <dt>Genomic Data Source</dt>
            <dd>VCF (version 4.2)</dd>
            <dt>In Model Version</dt>
            <dd>1</dd>
        </dl>
        <dl>
            <dt>Example</dt>
            <dd markdown="1">
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
</dd>
        </dl>
    </div>
</div>

<div class="panel panel-default">
    <div class="panel-heading">Key: effect</div>
    <div class="panel-body">
        <dl class="dl-horizontal">
            <dt>Description</dt>
            <dd>Name of the effect a variant has on another feature.</dd>
            <dt>Appears in</dt>
            <dd>feature object</dd>
            <dt>JSON Type</dt>
            <dd>string</dd>
            <dt>JSON-LD Type</dt>
            <dd>gfvo-squared:effect</dd>
            <dt>Genomic Data Source</dt>
            <dd>GVF (version 1.07)</dd>
            <dt>In Model Version</dt>
            <dd>1</dd>
        </dl>
        <dl>
            <dt>Example</dt>
            <dd markdown="1">
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
</dd>
        </dl>
    </div>
</div>

<div class="panel panel-default">
    <div class="panel-heading">Key: effects</div>
    <div class="panel-body">
        <dl class="dl-horizontal">
            <dt>Description</dt>
            <dd>Container for variant effects.</dd>
            <dt>Appears in</dt>
            <dd>feature object</dd>
            <dt>JSON Type</dt>
            <dd>array of objects</dd>
            <dt>JSON-LD Type</dt>
            <dd>gfvo-squared:effects</dd>
            <dt>Genomic Data Source</dt>
            <dd>GVF (version 1.07)</dd>
            <dt>In Model Version</dt>
            <dd>1</dd>
        </dl>
        <dl>
            <dt>Example</dt>
            <dd markdown="1">
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
</dd>
        </dl>
    </div>
</div>

<div class="panel panel-default">
    <div class="panel-heading">Key: end</div>
    <div class="panel-body">
        <dl class="dl-horizontal">
            <dt>Description</dt>
            <dd>End coordinate of a feature or alignment.</dd>
            <dt>Appears in</dt>
            <dd>feature object</dd>
            <dt>JSON Type</dt>
            <dd>number</dd>
            <dt>JSON-LD Type</dt>
            <dd>gfvo-squared:end</dd>
            <dt>Genomic Data Source</dt>
            <dd>GFF3 (version 1.21), GVF (version 1.07), VCF (version 4.2)</dd>
            <dt>In Model Version</dt>
            <dd>1</dd>
        </dl>
        <dl>
            <dt>Example</dt>
            <dd markdown="1">
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
</dd>
        </dl>
    </div>
</div>

<div class="panel panel-default">
    <div class="panel-heading">Key: experimentally-validated</div>
    <div class="panel-body">
        <dl class="dl-horizontal">
            <dt>Description</dt>
            <dd>Indicates whether sequence variant has been experimentally validated.</dd>
            <dt>Appears in</dt>
            <dd>feature object</dd>
            <dt>JSON Type</dt>
            <dd>boolean</dd>
            <dt>JSON-LD Type</dt>
            <dd>gfvo-squared:experimentallyValidated</dd>
            <dt>Genomic Data Source</dt>
            <dd>VCF (version 4.2)</dd>
            <dt>In Model Version</dt>
            <dd>1</dd>
        </dl>
        <dl>
            <dt>Example</dt>
            <dd markdown="1">
    {
        "@context" : "https://www.codamono.com/jsonld/vcf-f1.json",
        "experimentally-validated" : false
    }
</dd>
        </dl>
    </div>
</div>

<div class="panel panel-default">
    <div class="panel-heading">Key: feature-format</div>
    <div class="panel-body">
        <dl class="dl-horizontal">
            <dt>Description</dt>
            <dd>Information about key/values used to describe feature-centric data in a data set.</dd>
            <dt>Appears in</dt>
            <dd>meta object</dd>
            <dt>JSON Type</dt>
            <dd>object</dd>
            <dt>JSON-LD Type</dt>
            <dd>gfvo-squared:featureFormat</dd>
            <dt>Genomic Data Source</dt>
            <dd>VCF (version 4.2)</dd>
            <dt>In Model Version</dt>
            <dd>1</dd>
        </dl>
        <dl>
            <dt>Example</dt>
            <dd markdown="1">
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
</dd>
        </dl>
    </div>
</div>

<div class="panel panel-default">
    <div class="panel-heading">Key: features</div>
    <div class="panel-body">
        <dl class="dl-horizontal">
            <dt>Description</dt>
            <dd>Number of genomic-features lines that were read from a genomics file.</dd>
            <dt>Appears in</dt>
            <dd>summary object</dd>
            <dt>JSON Type</dt>
            <dd>number</dd>
            <dt>JSON-LD Type</dt>
            <dd>gfvo-squared:features</dd>
            <dt>In Model Version</dt>
            <dd>1</dd>
        </dl>
        <dl>
            <dt>Example:</dt>
            <dd markdown="1">
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
</dd>
        </dl>
    </div>
</div>

<div class="panel panel-default">
    <div class="panel-heading">Key: features-filtered</div>
    <div class="panel-body">
        <dl class="dl-horizontal">
            <dt>Description</dt>
            <dd>Number of genomic-features lines that were filtered through the Python API.</dd>
            <dt>Appears in</dt>
            <dd>summary object</dd>
            <dt>JSON Type</dt>
            <dd>number</dd>
            <dt>JSON-LD Type</dt>
            <dd>gfvo-squared:featuresFiltered</dd>
            <dt>In Model Version</dt>
            <dd>1</dd>
        </dl>
        <dl>
            <dt>Example:</dt>
            <dd markdown="1">
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
</dd>
        </dl>
    </div>
</div>

<div class="panel panel-default">
    <div class="panel-heading">Key: file-date</div>
    <div class="panel-body">
        <dl class="dl-horizontal">
            <dt>Description</dt>
            <dd>Creation date of the file that contains the represented genomics data.</dd>
            <dt>Appears in</dt>
            <dd>meta object</dd>
            <dt>JSON Type</dt>
            <dd>string</dd>
            <dt>JSON-LD Type</dt>
            <dd>gfvo-squared:fileDate</dd>
            <dt>Genomic Data Source</dt>
            <dd>GVF (version 1.07), VCF (version 4.2)</dd>
            <dt>In Model Version</dt>
            <dd>1</dd>
        </dl>
        <dl>
            <dt>Example</dt>
            <dd markdown="1">
    {
        "@context" : "https://www.codamono.com/jsonld/gvf-x1.json",
        "file-date" : "2015-03-08"
    }
</dd>
        </dl>
    </div>
</div>

<div class="panel panel-default">
    <div class="panel-heading">Key: finish</div>
    <div class="panel-body">
        <dl class="dl-horizontal">
            <dt>Description</dt>
            <dd>Finish time of the BioInterchange software.</dd>
            <dt>Appears in</dt>
            <dd>summary object</dd>
            <dt>JSON Type</dt>
            <dd>string</dd>
            <dt>JSON-LD Type</dt>
            <dd>gfvo-squared:finish</dd>
            <dt>In Model Version</dt>
            <dd>1</dd>
        </dl>
        <dl>
            <dt>Example:</dt>
            <dd markdown="1">
    {
        "@context" : "https://www.codamono.com/jsonld/biointerchange-s1.json",
         "runtime" : {
            "invocation" : "Fri Jul 17 17:00:54 2015",
            "finish" : "Fri Jul 17 17:00:56 2015",
            "lapsed-seconds" : "2"
        }
    }
</dd>
        </dl>
    </div>
</div>

<div class="panel panel-default">
    <div class="panel-heading">Key: genome-build</div>
    <div class="panel-body">
        <dl class="dl-horizontal">
            <dt>Description</dt>
            <dd>Information about the underlying genome builds of a data set.</dd>
            <dt>Appears in</dt>
            <dd>meta object</dd>
            <dt>JSON Type</dt>
            <dd>array of objects</dd>
            <dt>JSON-LD Type</dt>
            <dd>gfvo-squared:genomeBuild</dd>
            <dt>Genomic Data Source</dt>
            <dd>GFF3 (version 1.21), GVF (version 1.07)</dd>
            <dt>In Model Version</dt>
            <dd>1</dd>
        </dl>
        <dl>
            <dt>Example</dt>
            <dd markdown="1">
    {
        "@context" : "https://www.codamono.com/jsonld/gvf-f1.json",
        "genome-build" : [
            {
                "source" : "NCBI",
                "build" : "B36"
            }
        ]
    }
</dd>
        </dl>
    </div>
</div>

<div class="panel panel-default">
    <div class="panel-heading">Key: genomic-source</div>
    <div class="panel-body">
        <dl class="dl-horizontal">
            <dt>Description</dt>
            <dd>Information about the genomic origin source (LOINC code) of feature data.</dd>
            <dt>Appears in</dt>
            <dd>meta object</dd>
            <dt>JSON Type</dt>
            <dd>string</dd>
            <dt>JSON-LD Type</dt>
            <dd>gfvo-squared:genomicSource</dd>
            <dt>Genomic Data Source</dt>
            <dd>GVF (version 1.07)</dd>
            <dt>In Model Version</dt>
            <dd>1</dd>
        </dl>
        <dl>
            <dt>Example</dt>
            <dd markdown="1">
    {
        "@context" : "https://www.codamono.com/jsonld/gvf-x1.json",
        "genomic-source" : "somatic"
    }
</dd>
        </dl>
    </div>
</div>

<div class="panel panel-default">
    <div class="panel-heading">Key: genotype</div>
    <div class="panel-body">
        <dl class="dl-horizontal">
            <dt>Description</dt>
            <dd>Feature or sample specific genotype information.</dd>
            <dt>Appears in</dt>
            <dd>feature object</dd>
            <dt>JSON Type</dt>
            <dd>object</dd>
            <dt>JSON-LD Type</dt>
            <dd>gfvo-squared:genotype</dd>
            <dt>Genomic Data Source</dt>
            <dd>GVF (version 1.07), VCF (version 4.2)</dd>
            <dt>See Also</dt>
            <dd>AA, AB, etc.</dd>
            <dt>In Model Version</dt>
            <dd>1</dd>
        </dl>
        <dl>
            <dt>Example</dt>
            <dd markdown="1">
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
</dd>
        </dl>
    </div>
</div>

<div class="panel panel-default">
    <div class="panel-heading">Key: genotype-format</div>
    <div class="panel-body">
        <dl class="dl-horizontal">
            <dt>Description</dt>
            <dd>Information about keys/values used with genotypes in a data set.</dd>
            <dt>Appears in</dt>
            <dd>meta object</dd>
            <dt>JSON Type</dt>
            <dd>object</dd>
            <dt>JSON-LD Type</dt>
            <dd>gfvo-squared:genotypeFormat</dd>
            <dt>Genomic Data Source</dt>
            <dd>VCF (version 4.2)</dd>
            <dt>In Model Version</dt>
            <dd>1</dd>
        </dl>
        <dl>
            <dt>Example</dt>
            <dd markdown="1">
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
</dd>
        </dl>
    </div>
</div>

<div class="panel panel-default">
    <div class="panel-heading">Key: genotype-likelihood</div>
    <div class="panel-body">
        <dl class="dl-horizontal">
            <dt>Description</dt>
            <dd>Log-10 scaled genotype likelihood.</dd>
            <dt>Appears in</dt>
            <dd>feature object</dd>
            <dt>JSON Type</dt>
            <dd>number</dd>
            <dt>JSON-LD Type</dt>
            <dd>gfvo-squared:genotypeLikelihood</dd>
            <dt>Genomic Data Source</dt>
            <dd>VCF (version 4.2)</dd>
            <dt>See Also</dt>
            <dd>AA, AB, etc.</dd>
            <dt>In Model Version</dt>
            <dd>1</dd>
        </dl>
        <dl>
            <dt>Example</dt>
            <dd markdown="1">
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
</dd>
        </dl>
    </div>
</div>

<div class="panel panel-default">
    <div class="panel-heading">Key: genotype-likelihood-phred-scaled</div>
    <div class="panel-body">
        <dl class="dl-horizontal">
            <dt>Description</dt>
            <dd>Phred-scaled genotype likelihood.</dd>
            <dt>Appears in</dt>
            <dd>feature object</dd>
            <dt>JSON Type</dt>
            <dd>number</dd>
            <dt>JSON-LD Type</dt>
            <dd>gfvo-squared:genotypeLikelihoodPhredScaled</dd>
            <dt>Genomic Data Source</dt>
            <dd>VCF (version 4.2)</dd>
            <dt>See Also</dt>
            <dd>AA, AB, etc.</dd>
            <dt>In Model Version</dt>
            <dd>1</dd>
        </dl>
        <dl>
            <dt>Example</dt>
            <dd markdown="1">
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
</dd>
        </dl>
    </div>
</div>

<div class="panel panel-default">
    <div class="panel-heading">Key: gff-version</div>
    <div class="panel-body">
        <dl class="dl-horizontal">
            <dt>Description</dt>
            <dd>Genomics file-format versioning of GFF3 data sources.</dd>
            <dt>Appears in</dt>
            <dd>meta object</dd>
            <dt>JSON Type</dt>
            <dd>string</dd>
            <dt>JSON-LD Type</dt>
            <dd>gfvo-squared:gffVersion</dd>
            <dt>Genomic Data Source</dt>
            <dd>GFF3 (version 1.21), GVF (version 1.07)</dd>
            <dt>In Model Version</dt>
            <dd>1</dd>
        </dl>
        <dl>
            <dt>Example</dt>
            <dd markdown="1">
    {
        "@context" : "https://www.codamono.com/jsonld/gff3-x1.json",
        "gff-version" : "1.21"
    }
</dd>
        </dl>
    </div>
</div>

<div class="panel panel-default">
    <div class="panel-heading">Key: global</div>
    <div class="panel-body">
        <dl class="dl-horizontal">
            <dt>Description</dt>
            <dd>Global assignment that applies to all features in the data set.</dd>
            <dt>Appears in</dt>
            <dd>meta object</dd>
            <dt>JSON Type</dt>
            <dd>object</dd>
            <dt>JSON-LD Type</dt>
            <dd>gfvo-squared:global</dd>
            <dt>Genomic Data Source</dt>
            <dd>GFF3 (version 1.21), GVF (version 1.07), VCF (version 4.2)</dd>
            <dt>In Model Version</dt>
            <dd>1</dd>
        </dl>
        <dl>
            <dt>Example</dt>
            <dd markdown="1">
    {
        "@context" : "https://www.codamono.com/jsonld/gvf-x1.json",
        "technology-platform" : {
            "global" : {
                "comment" : "Preliminary data."
            }
        }
    }
</dd>
        </dl>
    </div>
</div>

<div class="panel panel-default">
    <div class="panel-heading">Key: gvf-version</div>
    <div class="panel-body">
        <dl class="dl-horizontal">
            <dt>Description</dt>
            <dd>Genomics file-format versioning of GVF data sources.</dd>
            <dt>Appears in</dt>
            <dd>meta object</dd>
            <dt>JSON Type</dt>
            <dd>string</dd>
            <dt>JSON-LD Type</dt>
            <dd>gfvo-squared:gvfVersion</dd>
            <dt>Genomic Data Source</dt>
            <dd>GVF (version 1.07)</dd>
            <dt>In Model Version</dt>
            <dd>1</dd>
        </dl>
        <dl>
            <dt>Example</dt>
            <dd markdown="1">
    {
        "@context" : "https://www.codamono.com/jsonld/gvf-x1.json",
        "gvf-version" : "1.07"
    }
</dd>
        </dl>
    </div>
</div>

<div class="panel panel-default">
    <div class="panel-heading">Key: id</div>
    <div class="panel-body">
        <dl class="dl-horizontal">
            <dt>Description</dt>
            <dd>An identifier.</dd>
            <dt>Appears in</dt>
            <dd>feature object, meta object</dd>
            <dt>JSON Type</dt>
            <dd>string</dd>
            <dt>JSON-LD Type</dt>
            <dd>gfvo-squared:id</dd>
            <dt>Genomic Data Source</dt>
            <dd>GFF3 (version 1.21), GVF (version 1.07), VCF (version 4.2)</dd>
            <dt>Note</dt>
            <dd>Multiple "ID" entries in VCF: only the first entry becomes an "id" in JSON-LD, where as the remaining identifiers become "alias" entries; this behaviour is in concordance with GFF3/GVF representations of identifiers.</dd>
            <dt>In Model Version</dt>
            <dd>1</dd>
        </dl>
        <dl>
            <dt>Example</dt>
            <dd markdown="1">
    {
        "@context" : "https://www.codamono.com/jsonld/gff3-f1.json",
        "id" : "ENSG00000139618"
    }
</dd>
        </dl>
    </div>
</div>

<div class="panel panel-default">
    <div class="panel-heading">Key: individuals</div>
    <div class="panel-body">
        <dl class="dl-horizontal">
            <dt>Description</dt>
            <dd>Identifiers of sequenced individuals.</dd>
            <dt>Appears in</dt>
            <dd>meta object</dd>
            <dt>JSON Type</dt>
            <dd>array of strings</dd>
            <dt>JSON-LD Type</dt>
            <dd>gfvo-squared:individuals</dd>
            <dt>Genomic Data Source</dt>
            <dd>GVF (version 1.07)</dd>
            <dt>In Model Version</dt>
            <dd>1</dd>
        </dl>
        <dl>
            <dt>Example</dt>
            <dd markdown="1">
    {
        "@context" : "https://www.codamono.com/jsonld/gvf-x1.json",
        "individuals" : [
            "NA18507",
            "NA12878",
            "NA19240"
        ]
    }
</dd>
        </dl>
    </div>
</div>

<div class="panel panel-default">
    <div class="panel-heading">Key: input-file</div>
    <div class="panel-body">
        <dl class="dl-horizontal">
            <dt>Description</dt>
            <dd>Name and relative path of the input file from which BioInterchange read data.</dd>
            <dt>Appears in</dt>
            <dd>context object</dd>
            <dt>JSON Type</dt>
            <dd>string</dd>
            <dt>JSON-LD Type</dt>
            <dd>gfvo-squared:inputFile</dd>
            <dt>In Model Version</dt>
            <dd>1</dd>
        </dl>
        <dl>
            <dt>Example</dt>
            <dd markdown="1">
    {
        "@context" : "https://www.codamono.com/jsonld/vcf-c1.json",
        "input-file" : "test-data/playground-vcf.vcf"
    }
</dd>
        </dl>
    </div>
</div>

<div class="panel panel-default">
    <div class="panel-heading">Key: input-filetype</div>
    <div class="panel-body">
        <dl class="dl-horizontal">
            <dt>Description</dt>
            <dd>File type of the input.</dd>
            <dt>Appears in</dt>
            <dd>context object</dd>
            <dt>JSON Type</dt>
            <dd>string (either &ldquo;GFF3&rdquo;, &ldquo;GVF&rdquo;, or &ldquo;VCF&rdquo;)</dd>
            <dt>JSON-LD Type</dt>
            <dd>gfvo-squared:inputFiletype</dd>
            <dt>In Model Version</dt>
            <dd>1</dd>
        </dl>
        <dl>
            <dt>Example</dt>
            <dd markdown="1">
    {
        "@context" : "https://www.codamono.com/jsonld/vcf-c1.json",
        "input-filetype" : "VCF"
    }
</dd>
        </dl>
    </div>
</div>

<div class="panel panel-default">
    <div class="panel-heading">Key: invocation</div>
    <div class="panel-body">
        <dl class="dl-horizontal">
            <dt>Description</dt>
            <dd>Invocation time of the BioInterchange software.</dd>
            <dt>Appears in</dt>
            <dd>summary object</dd>
            <dt>JSON Type</dt>
            <dd>string</dd>
            <dt>JSON-LD Type</dt>
            <dd>gfvo-squared:invocation</dd>
            <dt>In Model Version</dt>
            <dd>1</dd>
        </dl>
        <dl>
            <dt>Example:</dt>
            <dd markdown="1">
    {
        "@context" : "https://www.codamono.com/jsonld/biointerchange-s1.json",
         "runtime" : {
            "invocation" : "Fri Jul 17 17:00:54 2015",
            "finish" : "Fri Jul 17 17:00:56 2015",
            "lapsed-seconds" : "2"
        }
    }
</dd>
        </dl>
    </div>
</div>

<div class="panel panel-default">
    <div class="panel-heading">Key: landmark</div>
    <div class="panel-body">
        <dl class="dl-horizontal">
            <dt>Description</dt>
            <dd>A named genomic or proteomic landmark.</dd>
            <dt>Appears in</dt>
            <dd>feature object</dd>
            <dt>JSON Type</dt>
            <dd>string</dd>
            <dt>JSON-LD Type</dt>
            <dd>gfvo-squared:landmark</dd>
            <dt>Genomic Data Source</dt>
            <dd>GFF3 (version 1.21), GVF (version 1.07), VCF (version 4.2)</dd>
            <dt>In Model Version</dt>
            <dd>1</dd>
        </dl>
        <dl>
            <dt>Example</dt>
            <dd markdown="1">
    {
        "@context" : "https://www.codamono.com/jsonld/gff3-f1.json",
        "locus" : {
            "landmark" : "Chr1",
            "start" : 290207,
            "end" : 291002
        }
    }
</dd>
        </dl>
    </div>
</div>

<div class="panel panel-default">
    <div class="panel-heading">Key: landmarks</div>
    <div class="panel-body">
        <dl class="dl-horizontal">
            <dt>Description</dt>
            <dd>List of landmark identifiers.</dd>
            <dt>Appears in</dt>
            <dd>meta object</dd>
            <dt>JSON Type</dt>
            <dd>array of strings</dd>
            <dt>JSON-LD Type</dt>
            <dd>gfvo-squared:landmarks</dd>
            <dt>Genomic Data Source</dt>
            <dd>GVF (version 1.07)</dd>
            <dt>In Model Version</dt>
            <dd>1</dd>
        </dl>
        <dl>
            <dt>Example</dt>
            <dd markdown="1">
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
</dd>
        </dl>
    </div>
</div>

<div class="panel panel-default">
    <div class="panel-heading">Key: lapsed-seconds</div>
    <div class="panel-body">
        <dl class="dl-horizontal">
            <dt>Description</dt>
            <dd>Number of seconds that the BioInterchange software took for one execution.</dd>
            <dt>Appears in</dt>
            <dd>summary object</dd>
            <dt>JSON Type</dt>
            <dd>string</dd>
            <dt>JSON-LD Type</dt>
            <dd>gfvo-squared:lapsed-seconds</dd>
            <dt>In Model Version</dt>
            <dd>1</dd>
        </dl>
        <dl>
            <dt>Example:</dt>
            <dd markdown="1">
    {
        "@context" : "https://www.codamono.com/jsonld/biointerchange-s1.json",
         "runtime" : {
            "invocation" : "Fri Jul 17 17:00:54 2015",
            "finish" : "Fri Jul 17 17:00:56 2015",
            "lapsed-seconds" : "2"
        }
    }
</dd>
        </dl>
    </div>
</div>

<div class="panel panel-default">
    <div class="panel-heading">Key: length</div>
    <div class="panel-body">
        <dl class="dl-horizontal">
            <dt>Description</dt>
            <dd>Length of a genomic feature or continuous sequence region.</dd>
            <dt>Appears in</dt>
            <dd>meta object</dd>
            <dt>JSON Type</dt>
            <dd>number</dd>
            <dt>JSON-LD Type</dt>
            <dd>gfvo-squared:length</dd>
            <dt>Genomic Data Source</dt>
            <dd>GFF3 (version 1.21), GVF (version 1.07), VCF (version 4.2)</dd>
            <dt>In Model Version</dt>
            <dd>1</dd>
        </dl>
        <dl>
            <dt>Example</dt>
            <dd markdown="1">
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
</dd>
        </dl>
    </div>
</div>

<div class="panel panel-default">
    <div class="panel-heading">Key: locus</div>
    <div class="panel-body">
        <dl class="dl-horizontal">
            <dt>Description</dt>
            <dd>A genomic or proteomic locus.</dd>
            <dt>Appears in</dt>
            <dd>feature object</dd>
            <dt>JSON Type</dt>
            <dd>object</dd>
            <dt>JSON-LD Type</dt>
            <dd>gfvo-squared:locus</dd>
            <dt>Genomic Data Source</dt>
            <dd>GFF3 (version 1.21), GVF (version 1.07), VCF (version 4.2)</dd>
            <dt>In Model Version</dt>
            <dd>1</dd>
        </dl>
        <dl>
            <dt>Example</dt>
            <dd markdown="1">
    {
        "@context" : "https://www.codamono.com/jsonld/gff3-f1.json",
        "locus" : {
            "landmark" : "Chr1",
            "start" : 290207,
            "end" : 291002
        }
    }
</dd>
        </dl>
    </div>
</div>

<div class="panel panel-default">
    <div class="panel-heading">Key: mapping-quality-rms</div>
    <div class="panel-body">
        <dl class="dl-horizontal">
            <dt>Description</dt>
            <dd>Root-mean-square mapping quality.</dd>
            <dt>Appears in</dt>
            <dd>feature object</dd>
            <dt>JSON Type</dt>
            <dd>number</dd>
            <dt>JSON-LD Type</dt>
            <dd>gfvo-squared:mappingQualityRMS</dd>
            <dt>Genomic Data Source</dt>
            <dd>VCF (version 4.2)</dd>
            <dt>In Model Version</dt>
            <dd>1</dd>
        </dl>
        <dl>
            <dt>Example</dt>
            <dd markdown="1">
    {
        "@context" : "https://www.codamono.com/jsonld/vcf-f1.json",
        "mapping-quality-rms" : 29
    }
</dd>
        </dl>
    </div>
</div>

<div class="panel panel-default">
    <div class="panel-heading">Key: membership-1000G</div>
    <div class="panel-body">
        <dl class="dl-horizontal">
            <dt>Description</dt>
            <dd>Membership in the 1000 Genomes (1000G) project.</dd>
            <dt>Appears in</dt>
            <dd>feature object</dd>
            <dt>JSON Type</dt>
            <dd>boolean</dd>
            <dt>JSON-LD Type</dt>
            <dd>gfvo-squared:membership1000G</dd>
            <dt>Genomic Data Source</dt>
            <dd>VCF (version 4.2)</dd>
            <dt>In Model Version</dt>
            <dd>1</dd>
        </dl>
        <dl>
            <dt>Example</dt>
            <dd markdown="1">
    {
        "@context" : "https://www.codamono.com/jsonld/vcf-f1.json",
        "membership-1000G" : true
    }
</dd>
        </dl>
    </div>
</div>

<div class="panel panel-default">
    <div class="panel-heading">Key: membership-dbsnp</div>
    <div class="panel-body">
        <dl class="dl-horizontal">
            <dt>Description</dt>
            <dd>Membership in dbSNP.</dd>
            <dt>Appears in</dt>
            <dd>feature object</dd>
            <dt>JSON Type</dt>
            <dd>boolean</dd>
            <dt>JSON-LD Type</dt>
            <dd>gfvo-squared:membershipDbSNP</dd>
            <dt>Genomic Data Source</dt>
            <dd>VCF (version 4.2)</dd>
            <dt>In Model Version</dt>
            <dd>1</dd>
        </dl>
        <dl>
            <dt>Example</dt>
            <dd markdown="1">
    {
        "@context" : "https://www.codamono.com/jsonld/vcf-f1.json",
        "membership-dbsnp" : true
    }
</dd>
        </dl>
    </div>
</div>

<div class="panel panel-default">
    <div class="panel-heading">Key: membership-hapmap-2</div>
    <div class="panel-body">
        <dl class="dl-horizontal">
            <dt>Description</dt>
            <dd>Membership in HapMap 2.</dd>
            <dt>Appears in</dt>
            <dd>feature object</dd>
            <dt>JSON Type</dt>
            <dd>boolean</dd>
            <dt>JSON-LD Type</dt>
            <dd>gfvo-squared:membershipHapMap2</dd>
            <dt>Genomic Data Source</dt>
            <dd>VCF (version 4.2)</dd>
            <dt>In Model Version</dt>
            <dd>1</dd>
        </dl>
        <dl>
            <dt>Example</dt>
            <dd markdown="1">
    {
        "@context" : "https://www.codamono.com/jsonld/vcf-f1.json",
        "membership-hapmap-2" : true
    }
</dd>
        </dl>
    </div>
</div>

<div class="panel panel-default">
    <div class="panel-heading">Key: membership-hapmap-3</div>
    <div class="panel-body">
        <dl class="dl-horizontal">
            <dt>Description</dt>
            <dd>Membership in HapMap 3.
            <dt>Appears in</dt>
            <dd>feature object</dd>
            <dt>JSON Type</dt>
            <dd>boolean</dd>
            <dt>JSON-LD Type</dt>
            <dd>gfvo-squared:membershipHapMap3</dd>
            <dt>Genomic Data Source</dt>
            <dd>VCF (version 4.2)</dd>
            <dt>In Model Version</dt>
            <dd>1</dd>
        </dl>
        <dl>
            <dt>Example</dt>
            <dd markdown="1">
    {
        "@context" : "https://www.codamono.com/jsonld/vcf-f1.json",
        "membership-hapmap-3" : true
    }
</dd>
        </dl>
    </div>
</div>

<div class="panel panel-default">
    <div class="panel-heading">Key: meta-lines</div>
    <div class="panel-body">
        <dl class="dl-horizontal">
            <dt>Description</dt>
            <dd>Number of genomic-meta/pragma lines that were read from a genomics file.</dd>
            <dt>Appears in</dt>
            <dd>summary object</dd>
            <dt>JSON Type</dt>
            <dd>number</dd>
            <dt>JSON-LD Type</dt>
            <dd>gfvo-squared:metaLines</dd>
            <dt>In Model Version</dt>
            <dd>1</dd>
        </dl>
        <dl>
            <dt>Example:</dt>
            <dd markdown="1">
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
</dd>
        </dl>
    </div>
</div>

<div class="panel panel-default">
    <div class="panel-heading">Key: meta-lines-filtered</div>
    <div class="panel-body">
        <dl class="dl-horizontal">
            <dt>Description</dt>
            <dd>Number of genomic-features lines that were filtered through the Python API.</dd>
            <dt>Appears in</dt>
            <dd>summary object</dd>
            <dt>JSON Type</dt>
            <dd>number</dd>
            <dt>JSON-LD Type</dt>
            <dd>gfvo-squared:metaLinesFiltered</dd>
            <dt>In Model Version</dt>
            <dd>1</dd>
        </dl>
        <dl>
            <dt>Example:</dt>
            <dd markdown="1">
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
</dd>
        </dl>
    </div>
</div>

<div class="panel panel-default">
    <div class="panel-heading">Key: ontology</div>
    <div class="panel-body">
        <dl class="dl-horizontal">
            <dt>Description</dt>
            <dd>A reference to an external ontology.</dd>
            <dt>Appears in</dt>
            <dd>meta object</dd>
            <dt>JSON Type</dt>
            <dd>string</dd>
            <dt>JSON-LD Type</dt>
            <dd>gfvo-squared:ontology</dd>
            <dt>Genomic Data Source</dt>
            <dd>GFF3 (version 1.21), GVF (version 1.07)</dd>
            <dt>In Model Version</dt>
            <dd>1</dd>
        </dl>
        <dl>
            <dt>Example</dt>
            <dd markdown="1">
    {
        "@context" : "https://www.codamono.com/jsonld/gvf-x1.json",
        "phenotype-description" : {
            "global" : {
                "ontology" : "http://www.human-phenotype-ontology.org/human-phenotype-ontology.obo.gz"
            }
        }
    }
</dd>
        </dl>
    </div>
</div>

<div class="panel panel-default">
    <div class="panel-heading">Key: ontology-term</div>
    <div class="panel-body">
        <dl class="dl-horizontal">
            <dt>Description</dt>
            <dd>Ontology terms associated with a genomic feature.</dd>
            <dt>Appears in</dt>
            <dd>feature object</dd>
            <dt>JSON Type</dt>
            <dd>array of strings</dd>
            <dt>JSON-LD Type</dt>
            <dd>gfvo-squared:ontology-term</dd>
            <dt>Genomic Data Source</dt>
            <dd>GFF3 (version 1.21), GVF (version 1.07)</dd>
            <dt>In Model Version</dt>
            <dd>1</dd>
        </dl>
        <dl>
            <dt>Example</dt>
            <dd markdown="1">
    {
        "@context" : "https://www.codamono.com/jsonld/gff3-f1.json",
        "ontology-term" : [
            "GO:0046703"
        ]
    }
</dd>
        </dl>
    </div>
</div>

<div class="panel panel-default">
    <div class="panel-heading">Key: output-file</div>
    <div class="panel-body">
        <dl class="dl-horizontal">
            <dt>Description</dt>
            <dd>Filename and relative path to the output of the BioInterchange software.</dd>
            <dt>Appears in</dt>
            <dd>context object</dd>
            <dt>JSON Type</dt>
            <dd>string (&ldquo;null&rdquo;, if output was written to the console)</dd>
            <dt>JSON-LD Type</dt>
            <dd>gfvo-squared:outputFile</dd>
            <dt>In Model Version</dt>
            <dd>1</dd>
        </dl>
        <dl>
            <dt>Example</dt>
            <dd markdown="1">
    {
        "@context" : "https://www.codamono.com/jsonld/vcf-c1.json",
        "output-file" : "example.ldj"
    }
</dd>
        </dl>
    </div>
</div>

<div class="panel panel-default">
    <div class="panel-heading">Key: phase-set</div>
    <div class="panel-body">
        <dl class="dl-horizontal">
            <dt>Description</dt>
            <dd>Phase set identifier; indicates to which phase set a phased genotype belongs.</dd>
            <dt>Appears in</dt>
            <dd>feature object</dd>
            <dt>JSON Type</dt>
            <dd>number</dd>
            <dt>JSON-LD Type</dt>
            <dd>gfvo-squared:phaseSet</dd>
            <dt>Genomic Data Source</dt>
            <dd>VCF (version 4.2)</dd>
            <dt>In Model Version</dt>
            <dd>1</dd>
        </dl>
        <dl>
            <dt>Example</dt>
            <dd markdown="1">
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
</dd>
        </dl>
    </div>
</div>

<div class="panel panel-default">
    <div class="panel-heading">Key: phased</div>
    <div class="panel-body">
        <dl class="dl-horizontal">
            <dt>Description</dt>
            <dd>Indicates whether a genotype is phased or unphased.</dd>
            <dt>Appears in</dt>
            <dd>feature object</dd>
            <dt>JSON Type</dt>
            <dd>boolean</dd>
            <dt>JSON-LD Type</dt>
            <dd>gfvo-squared:phased</dd>
            <dt>Genomic Data Source</dt>
            <dd>GVF (version 1.07), VCF (version 4.2)</dd>
            <dt>In Model Version</dt>
            <dd>1</dd>
        </dl>
        <dl>
            <dt>Example</dt>
            <dd markdown="1">
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
</dd>
        </dl>
    </div>
</div>

<div class="panel panel-default">
    <div class="panel-heading">Key: phased-genotypes</div>
    <div class="panel-body">
        <dl class="dl-horizontal">
            <dt>Description</dt>
            <dd>Marks phased genotypes and provides extra information about them.</dd>
            <dt>Appears in</dt>
            <dd>meta object</dd>
            <dt>JSON Type</dt>
            <dd>boolean</dd>
            <dt>JSON-LD Type</dt>
            <dd>gfvo-squared:phasedGenotypes</dd>
            <dt>Genomic Data Source</dt>
            <dd>GVF (version 1.07)</dd>
            <dt>In Model Version</dt>
            <dd>1</dd>
        </dl>
        <dl>
            <dt>Example</dt>
            <dd markdown="1">
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
</dd>
        </dl>
    </div>
</div>

<div class="panel panel-default">
    <div class="panel-heading">Key: phasing-quality</div>
    <div class="panel-body">
        <dl class="dl-horizontal">
            <dt>Description</dt>
            <dd>Phred-scaled probability that alleles are ordered incorrectly.</dd>
            <dt>Appears in</dt>
            <dd>feature object</dd>
            <dt>JSON Type</dt>
            <dd>number</dd>
            <dt>JSON-LD Type</dt>
            <dd>gfvo-squared:phasingQuality</dd>
            <dt>Genomic Data Source</dt>
            <dd>VCF (version 4.2)</dd>
            <dt>In Model Version</dt>
            <dd>1</dd>
        </dl>
        <dl>
            <dt>Example</dt>
            <dd markdown="1">
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
</dd>
        </dl>
    </div>
</div>

<div class="panel panel-default">
    <div class="panel-heading">Key: phenotype-description</div>
    <div class="panel-body">
        <dl class="dl-horizontal">
            <dt>Description</dt>
            <dd>Additional information about phenotypes.</dd>
            <dt>Appears in</dt>
            <dd>meta object</dd>
            <dt>JSON Type</dt>
            <dd>object</dd>
            <dt>JSON-LD Type</dt>
            <dd>gfvo-squared:phenotypeDescription</dd>
            <dt>Genomic Data Source</dt>
            <dd>GVF (version 1.07)</dd>
            <dt>In Model Version</dt>
            <dd>1</dd>
        </dl>
        <dl>
            <dt>Example</dt>
            <dd markdown="1">
    {
        "@context" : "https://www.codamono.com/jsonld/gvf-x1.json",
        "phenotype-description" : {
            "global" : {
                "ontology" : "http://www.human-phenotype-ontology.org/human-phenotype-ontology.obo.gz"
            }
        }
    }
</dd>
        </dl>
    </div>
</div>

<div class="panel panel-default">
    <div class="panel-heading">Key: population</div>
    <div class="panel-body">
        <dl class="dl-horizontal">
            <dt>Description</dt>
            <dd>Population code that is assigned to an individual.</dd>
            <dt>Appears in</dt>
            <dd>meta object</dd>
            <dt>JSON Type</dt>
            <dd>string</dd>
            <dt>JSON-LD Type</dt>
            <dd>gfvo-squared:population</dd>
            <dt>Genomic Data Source</dt>
            <dd>GVF (version 1.07)</dd>
            <dt>See Also</dt>
            <dd><a href="ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/README.populations">1000 Genomes project population codes</a></dd>
            <dt>In Model Version</dt>
            <dd>1</dd>
        </dl>
        <dl>
            <dt>Example</dt>
            <dd markdown="1">
    {
        "@context" : "https://www.codamono.com/jsonld/gvf-x1.json",
        "population" : "YRI"
    }
</dd>
        </dl>
    </div>
</div>

<div class="panel panel-default">
    <div class="panel-heading">Key: python-callback</div>
    <div class="panel-body">
        <dl class="dl-horizontal">
            <dt>Description</dt>
            <dd>Package and module name used for the Python API.</dd>
            <dt>Appears in</dt>
            <dd>context object</dd>
            <dt>JSON Type</dt>
            <dd>string (&ldquo;null&rdquo;, if the Python API was not used)</dd>
            <dt>JSON-LD Type</dt>
            <dd>gfvo-squared:pythonCallback</dd>
            <dt>In Model Version</dt>
            <dd>1</dd>
        </dl>
        <dl>
            <dt>Example:</dt>
            <dd markdown="1">
    {
        "@context" : "https://www.codamono.com/jsonld/vcf-c1.json",
        "python-callback" : "simplepy.simple"
    }
</dd>
        </dl>
    </div>
</div>

<div class="panel panel-default">
    <div class="panel-heading">Key: read-pair-span</div>
    <div class="panel-body">
        <dl class="dl-horizontal">
            <dt>Description</dt>
            <dd>Global assignment that applies to all features in the data set.</dd>
            <dt>Appears in</dt>
            <dd>meta object</dd>
            <dt>JSON Type</dt>
            <dd>array of numbers</dd>
            <dt>JSON-LD Type</dt>
            <dd>gfvo-squared:readPairSpan</dd>
            <dt>Genomic Data Source</dt>
            <dd>GVF (version 1.07)</dd>
            <dt>In Model Version</dt>
            <dd>1</dd>
        </dl>
        <dl>
            <dt>Example</dt>
            <dd markdown="1">
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
</dd>
        </dl>
    </div>
</div>

<div class="panel panel-default">
    <div class="panel-heading">Key: reads-with-zero-mapping-quality</div>
    <div class="panel-body">
        <dl class="dl-horizontal">
            <dt>Description</dt>
            <dd>Number of reads with mapping quality being equal to zero.</dd>
            <dt>Appears in</dt>
            <dd>feature object</dd>
            <dt>JSON Type</dt>
            <dd>number</dd>
            <dt>JSON-LD Type</dt>
            <dd>gfvo-squared:readsWithZeroMappingQuality</dd>
            <dt>Genomic Data Source</dt>
            <dd>VCF (version 4.2)</dd>
            <dt>In Model Version</dt>
            <dd>1</dd>
        </dl>
        <dl>
            <dt>Example</dt>
            <dd markdown="1">
    {
        "@context" : "https://www.codamono.com/jsonld/vcf-f1.json",
        "reads-with-zero-mapping-quality" : 2
    }
</dd>
        </dl>
    </div>
</div>

<div class="panel panel-default">
    <div class="panel-heading">Key: reference</div>
    <div class="panel-body">
        <dl class="dl-horizontal">
            <dt>Description</dt>
            <dd>Reference sequence information.</dd>
            <dt>Appears in</dt>
            <dd>feature object</dd>
            <dt>JSON Type</dt>
            <dd>object</dd>
            <dt>JSON-LD Type</dt>
            <dd>gfvo-squared:reference</dd>
            <dt>Genomic Data Source</dt>
            <dd>GVF (version 1.07), VCF (version 4.2)</dd>
            <dt>In Model Version</dt>
            <dd>1</dd>
        </dl>
        <dl>
            <dt>Example</dt>
            <dd markdown="1">
    {
        "@context" : "https://www.codamono.com/jsonld/vcf-f1.json",
        "reference" : {
            "sequence" : "A"
        }
    }
</dd>
        </dl>
    </div>
</div>

<div class="panel panel-default">
    <div class="panel-heading">Key: reference-fasta</div>
    <div class="panel-body">
        <dl class="dl-horizontal">
            <dt>Description</dt>
            <dd>A link to a FASTA file that contains reference sequences.</dd>
            <dt>Appears in</dt>
            <dd>meta object</dd>
            <dt>JSON Type</dt>
            <dd>string</dd>
            <dt>JSON-LD Type</dt>
            <dd>gfvo-squared:referenceFASTA</dd>
            <dt>Genomic Data Source</dt>
            <dd>GVF (version 1.07), VCF (version 4.2)</dd>
            <dt>In Model Version</dt>
            <dd>1</dd>
        </dl>
        <dl>
            <dt>Example</dt>
            <dd markdown="1">
    {
        "@context" : "https://www.codamono.com/jsonld/vcf-x1.json",
        "reference-fasta" : "ftp://ftp-mouse.sanger.ac.uk/ref/GRCm38_68.fa"
    }
</dd>
        </dl>
    </div>
</div>

<div class="panel panel-default">
    <div class="panel-heading">Key: runtime</div>
    <div class="panel-body">
        <dl class="dl-horizontal">
            <dt>Description</dt>
            <dd>Information about runtime statistics of the BioInterchange software.</dd>
            <dt>Appears in</dt>
            <dd>summary object</dd>
            <dt>JSON Type</dt>
            <dd>object</dd>
            <dt>JSON-LD Type</dt>
            <dd>gfvo-squared:runtime</dd>
            <dt>In Model Version</dt>
            <dd>1</dd>
        </dl>
        <dl>
            <dt>Example:</dt>
            <dd markdown="1">
    {
        "@context" : "https://www.codamono.com/jsonld/biointerchange-s1.json",
         "runtime" : {
            "invocation" : "Fri Jul 17 17:00:54 2015",
            "finish" : "Fri Jul 17 17:00:56 2015",
            "lapsed-seconds" : "2"
        }
    }
</dd>
        </dl>
    </div>
</div>

<div class="panel panel-default">
    <div class="panel-heading">Key: samples</div>
    <div class="panel-body">
        <dl class="dl-horizontal">
            <dt>Description</dt>
            <dd>Sample information</dd>
            <dt>Appears in</dt>
            <dd>feature object</dd>
            <dt>JSON Type</dt>
            <dd>orray of bject</dd>
            <dt>JSON-LD Type</dt>
            <dd>gfvo-squared:samples</dd>
            <dt>Genomic Data Source</dt>
            <dd>VCF (version 4.2)</dd>
            <dt>In Model Version</dt>
            <dd>1</dd>
        </dl>
        <dl>
            <dt>Example</dt>
            <dd markdown="1">
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
</dd>
        </dl>
    </div>
</div>

<div class="panel panel-default">
    <div class="panel-heading">Key: samples-with-data</div>
    <div class="panel-body">
        <dl class="dl-horizontal">
            <dt>Description</dt>
            <dd>Number of samples with data.</dd>
            <dt>Appears in</dt>
            <dd>feature object</dd>
            <dt>JSON Type</dt>
            <dd>number</dd>
            <dt>JSON-LD Type</dt>
            <dd>gfvo-squared:samplesWithData</dd>
            <dt>Genomic Data Source</dt>
            <dd>VCF (version 4.2)</dd>
            <dt>In Model Version</dt>
            <dd>1</dd>
        </dl>
        <dl>
            <dt>Example</dt>
            <dd markdown="1">
    {
        "@context" : "https://www.codamono.com/jsonld/vcf-f1.json",
        "samples-with-data" : 1
    }
</dd>
        </dl>
    </div>
</div>

<div class="panel panel-default">
    <div class="panel-heading">Key: score-method</div>
    <div class="panel-body">
        <dl class="dl-horizontal">
            <dt>Description</dt>
            <dd>Information about scoring methods that are used in a data set.</dd>
            <dt>Appears in</dt>
            <dd>meta object</dd>
            <dt>JSON Type</dt>
            <dd>object</dd>
            <dt>JSON-LD Type</dt>
            <dd>gfvo-squared:scoreMethod</dd>
            <dt>Genomic Data Source</dt>
            <dd>GVF (version 1.07)</dd>
            <dt>In Model Version</dt>
            <dd>1</dd>
        </dl>
        <dl>
            <dt>Example</dt>
            <dd markdown="1">
    {
        "@context" : "https://www.codamono.com/jsonld/gvf-x1.json",
        "score-method" : {
            "global" : {
                "comment" : "Scores are Phred scaled probabilities of an incorrect sequence_alteration call"
            }
        }
    }
</dd>
        </dl>
    </div>
</div>

<div class="panel panel-default">
    <div class="panel-heading">Key: sequence</div>
    <div class="panel-body">
        <dl class="dl-horizontal">
            <dt>Description</dt>
            <dd>Genomic or proteomic sequence.</dd>
            <dt>Appears in</dt>
            <dd>feature object</dd>
            <dt>JSON Type</dt>
            <dd>string</dd>
            <dt>JSON-LD Type</dt>
            <dd>gfvo-squared:sequence</dd>
            <dt>Genomic Data Source</dt>
            <dd>GFF3 (version 1.21), GVF (version 1.07), VCF (version 4.2)</dd>
            <dt>In Model Version</dt>
            <dd>1</dd>
        </dl>
        <dl>
            <dt>Example</dt>
            <dd markdown="1">
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
</dd>
        </dl>
    </div>
</div>

<div class="panel panel-default">
    <div class="panel-heading">Key: sequences</div>
    <div class="panel-body">
        <dl class="dl-horizontal">
            <dt>Description</dt>
            <dd>A list of genomic or proteomic sequences.</dd>
            <dt>Appears in</dt>
            <dd>feature object</dd>
            <dt>JSON Type</dt>
            <dd>object</dd>
            <dt>JSON-LD Type</dt>
            <dd>gfvo-squared:sequences</dd>
            <dt>Genomic Data Source</dt>
            <dd>VCF (version 4.2)</dd>
            <dt>In Model Version</dt>
            <dd>1</dd>
        </dl>
        <dl>
            <dt>Example</dt>
            <dd markdown="1">
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
</dd>
        </dl>
    </div>
</div>

<div class="panel panel-default">
    <div class="panel-heading">Key: sex</div>
    <div class="panel-body">
        <dl class="dl-horizontal">
            <dt>Description</dt>
            <dd>Sex of a sequenced individual.</dd>
            <dt>Appears in</dt>
            <dd>meta object</dd>
            <dt>JSON Type</dt>
            <dd>string</dd>
            <dt>JSON-LD Type</dt>
            <dd>gfvo-squared:sex</dd>
            <dt>Genomic Data Source</dt>
            <dd>GVF (version 1.07)</dd>
            <dt>In Model Version</dt>
            <dd>1</dd>
        </dl>
        <dl>
            <dt>Example</dt>
            <dd markdown="1">
    {
        "@context" : "https://www.codamono.com/jsonld/gvf-x1.json",
        "sex" : "male"
    }
</dd>
        </dl>
    </div>
</div>

<div class="panel panel-default">
    <div class="panel-heading">Key: somatic-mutation</div>
    <div class="panel-body">
        <dl class="dl-horizontal">
            <dt>Description</dt>
            <dd>Indicates whether a feature is a somatic mutation.</dd>
            <dt>Appears in</dt>
            <dd>feature object</dd>
            <dt>JSON Type</dt>
            <dd>boolean</dd>
            <dt>JSON-LD Type</dt>
            <dd>gfvo-squared:somaticMutation</dd>
            <dt>Genomic Data Source</dt>
            <dd>VCF (version 4.2)</dd>
            <dt>In Model Version</dt>
            <dd>1</dd>
        </dl>
        <dl>
            <dt>Example</dt>
            <dd markdown="1">
    {
        "@context" : "https://www.codamono.com/jsonld/vcf-f1.json",
        "somatic-mutation" : true
    }
</dd>
        </dl>
    </div>
</div>

<div class="panel panel-default">
    <div class="panel-heading">Key: source</div>
    <div class="panel-body">
        <dl class="dl-horizontal">
            <dt>Description</dt>
            <dd>Database or algorithm name that is the source of a feature or data set.</dd>
            <dt>Appears in</dt>
            <dd>feature object</dd>
            <dt>JSON Type</dt>
            <dd>string</dd>
            <dt>JSON-LD Type</dt>
            <dd>gfvo-squared:source</dd>
            <dt>Genomic Data Source</dt>
            <dd>GFF3 (version 1.21), GVF (version 1.07)</dd>
            <dt>In Model Version</dt>
            <dd>1</dd>
        </dl>
        <dl>
            <dt>Example</dt>
            <dd markdown="1">
    {
        "@context" : "https://www.codamono.com/jsonld/gvf-f1.json",
        "source" : "samtools"
    }
</dd>
        </dl>
    </div>
</div>

<div class="panel panel-default">
    <div class="panel-heading">Key: source-method</div>
    <div class="panel-body">
        <dl class="dl-horizontal">
            <dt>Description</dt>
            <dd>Information about data sources in a data set.</dd>
            <dt>Appears in</dt>
            <dd>meta object</dd>
            <dt>JSON Type</dt>
            <dd>object</dd>
            <dt>JSON-LD Type</dt>
            <dd>gfvo-squared:sourceMethod</dd>
            <dt>Genomic Data Source</dt>
            <dd>GVF (version 1.07)</dd>
            <dt>In Model Version</dt>
            <dd>1</dd>
        </dl>
        <dl>
            <dt>Example</dt>
            <dd markdown="1">
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
</dd>
        </dl>
    </div>
</div>

<div class="panel panel-default">
    <div class="panel-heading">Key: sources</div>
    <div class="panel-body">
        <dl class="dl-horizontal">
            <dt>Description</dt>
            <dd>List of data sources.</dd>
            <dt>Appears in</dt>
            <dd>meta object</dd>
            <dt>JSON Type</dt>
            <dd>array of strings</dd>
            <dt>JSON-LD Type</dt>
            <dd>gfvo-squared:sources</dd>
            <dt>Genomic Data Source</dt>
            <dd>GVF (version 1.07)</dd>
            <dt>In Model Version</dt>
            <dd>1</dd>
        </dl>
        <dl>
            <dt>Example</dt>
            <dd markdown="1">
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
</dd>
        </dl>
    </div>
</div>

<div class="panel panel-default">
    <div class="panel-heading">Key: start</div>
    <div class="panel-body">
        <dl class="dl-horizontal">
            <dt>Description</dt>
            <dd>Start coordinate of a feature or alignment.</dd>
            <dt>Appears in</dt>
            <dd>feature object</dd>
            <dt>JSON Type</dt>
            <dd>number</dd>
            <dt>JSON-LD Type</dt>
            <dd>gfvo-squared:start</dd>
            <dt>Genomic Data Source</dt>
            <dd>GFF3 (version 1.21), GVF (version 1.07), VCF (version 4.2)</dd>
            <dt>In Model Version</dt>
            <dd>1</dd>
        </dl>
        <dl>
            <dt>Example</dt>
            <dd markdown="1">
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
</dd>
        </dl>
    </div>
</div>

<div class="panel panel-default">
    <div class="panel-heading">Key: statistics</div>
    <div class="panel-body">
        <dl class="dl-horizontal">
            <dt>Description</dt>
            <dd>Summary statistics about the genomic data that was processed.</dd>
            <dt>Appears in</dt>
            <dd>summary object</dd>
            <dt>JSON Type</dt>
            <dd>object</dd>
            <dt>JSON-LD Type</dt>
            <dd>gfvo-squared:statistics</dd>
            <dt>In Model Version</dt>
            <dd>1</dd>
        </dl>
        <dl>
            <dt>Example:</dt>
            <dd markdown="1">
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
</dd>
        </dl>
    </div>
</div>

<div class="panel panel-default">
    <div class="panel-heading">Key: strand</div>
    <div class="panel-body">
        <dl class="dl-horizontal">
            <dt>Description</dt>
            <dd>Strand on which a feature is located.</dd>
            <dt>Appears in</dt>
            <dd>feature object</dd>
            <dt>JSON Type</dt>
            <dd>string ("+": forward strand, "-": reverse strand)</dd>
            <dt>JSON-LD Type</dt>
            <dd>gfvo-squared:strand</dd>
            <dt>Genomic Data Source</dt>
            <dd>GFF3 (version 1.21), GVF (version 1.07), VCF (version 4.2)</dd>
            <dt>In Model Version</dt>
            <dd>1</dd>
        </dl>
        <dl>
            <dt>Example</dt>
            <dd markdown="1">
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
</dd>
        </dl>
    </div>
</div>

<div class="panel panel-default">
    <div class="panel-heading">Key: strand-bias</div>
    <div class="panel-body">
        <dl class="dl-horizontal">
            <dt>Description</dt>
            <dd>Strand bias.</dd>
            <dt>Appears in</dt>
            <dd>feature object</dd>
            <dt>JSON Type</dt>
            <dd>number</dd>
            <dt>JSON-LD Type</dt>
            <dd>gfvo-squared:strandBias</dd>
            <dt>Genomic Data Source</dt>
            <dd>VCF (version 4.2)</dd>
            <dt>In Model Version</dt>
            <dd>1</dd>
        </dl>
        <dl>
            <dt>Example</dt>
            <dd markdown="1">
    {
        "@context" : "https://www.codamono.com/jsonld/vcf-f1.json",
        "strand-bias" : 0.5
    }
</dd>
        </dl>
    </div>
</div>

<div class="panel panel-default">
    <div class="panel-heading">Key: technology-platform</div>
    <div class="panel-body">
        <dl class="dl-horizontal">
            <dt>Description</dt>
            <dd>Information about the technology platform that was used to create a data set.</dd>
            <dt>Appears in</dt>
            <dd>meta object</dd>
            <dt>JSON Type</dt>
            <dd>object</dd>
            <dt>JSON-LD Type</dt>
            <dd>gfvo-squared:technologyPlatform</dd>
            <dt>Genomic Data Source</dt>
            <dd>GVF (version 1.07)</dd>
            <dt>In Model Version</dt>
            <dd>1</dd>
        </dl>
        <dl>
            <dt>Example</dt>
            <dd markdown="1">
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
</dd>
        </dl>
    </div>
</div>

<div class="panel panel-default">
    <div class="panel-heading">Key: type</div>
    <div class="panel-body">
        <dl class="dl-horizontal">
            <dt>Description</dt>
            <dd>Feature type.</dd>
            <dt>Appears in</dt>
            <dd>feature object</dd>
            <dt>JSON Type</dt>
            <dd>string</dd>
            <dt>JSON-LD Type</dt>
            <dd>gfvo-squared:type</dd>
            <dt>Genomic Data Source</dt>
            <dd>GFF3 (version 1.21), GVF (version 1.07)</dd>
            <dt>In Model Version</dt>
            <dd>1</dd>
        </dl>
        <dl>
            <dt>Example</dt>
            <dd markdown="1">
    {
        "@context" : "https://www.codamono.com/jsonld/gvf-f1.json",
        "type" : "SNV"
    }
</dd>
        </dl>
    </div>
</div>

<div class="panel panel-default">
    <div class="panel-heading">Key: types</div>
    <div class="panel-body">
        <dl class="dl-horizontal">
            <dt>Description</dt>
            <dd>List of feature types.</dd>
            <dt>Appears in</dt>
            <dd>meta object</dd>
            <dt>JSON Type</dt>
            <dd>array of strings</dd>
            <dt>JSON-LD Type</dt>
            <dd>gfvo-squared:types</dd>
            <dt>Genomic Data Source</dt>
            <dd>GVF (version 1.07)</dd>
            <dt>In Model Version</dt>
            <dd>1</dd>
        </dl>
        <dl>
            <dt>Example</dt>
            <dd markdown="1">
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
</dd>
        </dl>
    </div>
</div>

<div class="panel panel-default">
    <div class="panel-heading">Key: user-defined</div>
    <div class="panel-body">
        <dl class="dl-horizontal">
            <dt>Description</dt>
            <dd>Key/value pairs that are not defined in the GFF3-, GVF-, or VCF-specification.</dd>
            <dt>Appears in</dt>
            <dd>meta object, feature object</dd>
            <dt>JSON Type</dt>
            <dd>object containing key/value pairs</dd>
            <dt>JSON-LD Type</dt>
            <dd>gfvo-squared:user-defined</dd>
            <dt>Genomic Data Source</dt>
            <dd>GFF3 (version 1.21), GVF (version 1.07), VCF (version 4.2)</dd>
            <dt>See Also</dt>
            <dd>Appendix: unknown keys</dd>
            <dt>In Model Version</dt>
            <dd>1</dd>
        </dl>
        <dl>
            <dt>Example:</dt>
            <dd markdown="1">
    {
        "@context" : "https://www.codamono.com/jsonld/vcf-f1.json",
        "user-defined" : {
            "VDB" : 0.0006,
            "INDEL" : true,
            "DP4" : "0,0,71,0"
        }
    }
</dd>
        </dl>
    </div>
</div>

<div class="panel panel-default">
    <div class="panel-heading">Key: user-parameter</div>
    <div class="panel-body">
        <dl class="dl-horizontal">
            <dt>Description</dt>
            <dd>&ldquo;-u&rdquo; parameter value when running the BioInterchange software.</dd>
            <dt>Appears in</dt>
            <dd>context object</dd>
            <dt>JSON Type</dt>
            <dd>string (&ldquo;null&rdquo; when no &ldquo;-u&rdquo; was used)</dd>
            <dt>JSON-LD Type</dt>
            <dd>gfvo-squared:user-parameter</dd>
            <dt>See Also</dt>
            <dd>Appendix: unknown keys</dd>
            <dt>In Model Version</dt>
            <dd>1</dd>
        </dl>
        <dl>
            <dt>Example:</dt>
            <dd markdown="1">
    {
        "@context" : "https://www.codamono.com/jsonld/gff3-c1.json",
        "user-parameter" : "Preliminary data."
    }
</dd>
        </dl>
    </div>
</div>

<div class="panel panel-default">
    <div class="panel-heading">Key: variants</div>
    <div class="panel-body">
        <dl class="dl-horizontal">
            <dt>Description</dt>
            <dd>A collection of sequence variants.</dd>
            <dt>Appears in</dt>
            <dd>feature object</dd>
            <dt>JSON Type</dt>
            <dd>object</dd>
            <dt>JSON-LD Type</dt>
            <dd>gfvo-squared:variants</dd>
            <dt>Genomic Data Source</dt>
            <dd>GVF (version 1.07), VCF (version 4.2)</dd>
            <dt>In Model Version</dt>
            <dd>1</dd>
        </dl>
        <dl>
            <dt>Example</dt>
            <dd markdown="1">
    {
        "@context" : "https://www.codamono.com/jsonld/gvf-f1.json",
        "variants" : {
            "B" : {
                "sequence" : "T"
            }
        }
    }
</dd>
        </dl>
    </div>
</div>

<div class="panel panel-default">
    <div class="panel-heading">Key: vcf-version</div>
    <div class="panel-body">
        <dl class="dl-horizontal">
            <dt>Description</dt>
            <dd>Genomics file-format versioning of VCF data sources.</dd>
            <dt>Appears in</dt>
            <dd>meta object</dd>
            <dt>JSON Type</dt>
            <dd>string</dd>
            <dt>JSON-LD Type</dt>
            <dd>gfvo-squared:vcfVersion</dd>
            <dt>Genomic Data Source</dt>
            <dd>VCF (version 4.2)</dd>
            <dt>In Model Version</dt>
            <dd>1</dd>
        </dl>
        <dl>
            <dt>Example</dt>
            <dd markdown="1">
    {
        "@context" : "https://www.codamono.com/jsonld/vcf-x1.json",
        "vcf-version" : "4.2"
    }
</dd>
        </dl>
    </div>
</div>

<div class="panel panel-default">
    <div class="panel-heading">Appendix: unknown keys</div>
    <div class="panel-body">
        <dl class="dl-horizontal">
            <dt>Description</dt>
            <dd>Keys that are user defined.</dd>
            <dt>Appears in</dt>
            <dd>feature object, meta object</dd>
            <dt>JSON Type</dt>
            <dd>string</dd>
            <dt>JSON-LD Type</dt>
            <dd>gfvo-squared:unknownProperty-*</dd>
            <dt>Genomic Data Source</dt>
            <dd>GFF3 (version 1.21), GVF (version 1.07), VCF (version 4.2)</dd>
            <dt>In Model Version</dt>
            <dd>1</dd>
            <dt>Note</dt>
            <dd>The GFVO<sup>2</sup> type expands to include the unknown key. For example, "VDB" becomes "gfvo-squared:unknownProperty-VDB".</dd>
        </dl>
        <dl>
            <dt>Example</dt>
            <dd markdown="1">
    {
        "@context" : "https://www.codamono.com/jsonld/vcf-f1.json",
        "user-defined" : {
            "VDB" : 0.0006,
            "INDEL" : true,
            "DP4" : "0,0,71,0"
        }
    }
</dd>
        </dl>
    </div>
</div>

</div>

