"""
Extension for sourmash.sourmash_args
"""
import sys
import os
from sourmash.logging import notify, error
from sourmash.index import LinearIndex
from sourmash import signature as sig
from sourmash.sbt import SBT
from sourmash.sbtmh import SigLeaf
from sourmash.lca import LCA_Database
from sourmash.sourmash_args import get_moltype, traverse_find_sigs, filter_compatible_signatures, \
    check_tree_is_compatible

DEFAULT_LOAD_K = 31


def load_dbs_and_sigs(filenames, query, is_similarity_query, traverse=False):
    """
    Load one or more SBTs, LCAs, and/or signatures.

    Check for compatibility with query.
    """
    query_ksize = query.minhash.ksize
    query_moltype = get_moltype(query)

    n_signatures = 0
    n_databases = 0
    databases = []
    for sbt_or_sigfile in filenames:
        notify('loading from {}...', sbt_or_sigfile, end='\r')
        # are we collecting signatures from a directory/path?
        if traverse and os.path.isdir(sbt_or_sigfile):
            for sigfile in traverse_find_sigs([sbt_or_sigfile]):
                try:
                    siglist = sig.load_signatures(sigfile,
                                                  ksize=query_ksize,
                                                  select_moltype=query_moltype)
                    siglist = filter_compatible_signatures(query, siglist, 1)
                    linear = LinearIndex(siglist, filename=sigfile)
                    databases.append((linear, sbt_or_sigfile, False))
                    notify('loaded {} signatures from {}', len(linear),
                           sigfile, end='\r')
                    n_signatures += len(linear)
                except Exception:  # ignore errors with traverse
                    pass

            # done! jump to beginning of main 'for' loop
            continue

        # no traverse? try loading as an SBT.
        try:
            tree = SBT.load(sbt_or_sigfile, leaf_loader=SigLeaf.load)

            if not check_tree_is_compatible(sbt_or_sigfile, tree, query,
                                            is_similarity_query):
                sys.exit(-1)

            databases.append((tree, sbt_or_sigfile, 'SBT'))
            notify('loaded SBT {}', sbt_or_sigfile, end='\r')
            n_databases += 1

            # done! jump to beginning of main 'for' loop
            continue
        except (ValueError, EnvironmentError):
            # not an SBT - try as an LCA
            pass

        # ok. try loading as an LCA.
        try:
            lca_db = LCA_Database.load(sbt_or_sigfile)

            assert query_ksize == lca_db.ksize
            query_scaled = query.minhash.scaled

            notify('loaded LCA {}', sbt_or_sigfile, end='\r')
            n_databases += 1

            databases.append((lca_db, sbt_or_sigfile, 'LCA'))

            continue
        except (ValueError, TypeError, EnvironmentError):
            # not an LCA database - try as a .sig
            pass

        # not a tree? try loading as a signature.
        try:
            siglist = sig.load_signatures(sbt_or_sigfile,
                                          ksize=query_ksize,
                                          select_moltype=query_moltype)
            siglist = list(siglist)
            if len(siglist) == 0:  # file not found, or parse error?
                raise ValueError

            siglist = filter_compatible_signatures(query, siglist, False)
            linear = LinearIndex(siglist, filename=sbt_or_sigfile)
            databases.append((linear, sbt_or_sigfile, 'signature'))

            notify('loaded {} signatures from {}', len(linear),
                   sbt_or_sigfile, end='\r')
            n_signatures += len(linear)
        except (EnvironmentError, ValueError) as e:
            raise e

    if n_signatures and n_databases:
        notify('loaded {} signatures and {} databases total.', n_signatures,
               n_databases)
    elif n_signatures:
        notify('loaded {} signatures.', n_signatures)
    elif n_databases:
        notify('loaded {} databases.', n_databases)
    else:
        sys.exit(-1)

    if databases:
        print('')

    return databases
