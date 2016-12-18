#ifndef SURF_CONFIG_HPP
#define SURF_CONFIG_HPP

#include "sdsl/config.hpp"
#include <string>
#include <vector>

namespace surf {

const std::string TEXT_FILENAME = "text_int_SURF.sdsl";
const std::string TEXT_FILENAME_BYTE = "text_SURF.sdsl";
const std::string DICT_FILENAME = "dict.txt";
const std::string URL2ID_FILENAME = "url2id.txt";
const std::string DOCNAMES_FILENAME = "doc_names.txt";
const std::string SPACEUSAGE_FILENAME = "space_usage";

const std::string KEY_DOCWEIGHT = "docweights";
const std::string KEY_DARRAY = "darray";
const std::string KEY_U = "U";
const std::string KEY_WTU = "wtu";
const std::string KEY_DOCPERM = "docperm";
const std::string KEY_SADADF = "sadadf";
const std::string KEY_SADADF_G = "sadadf_g";
const std::string KEY_WTD = "wtd";
const std::string KEY_C = "C";
const std::string KEY_RMQC = "RMQC";
const std::string KEY_RMQW = "RMQW";
const std::string KEY_WTC = "wtc";
const std::string KEY_TMPCST = "tempcst";
const std::string KEY_TMPDUP = "tmpdup";
const std::string KEY_WTDUP  = "wtdup";
const std::string KEY_WTDP  = "wtdp";
const std::string KEY_DOC_OFFSET  = "doc_offset";
const std::string KEY_DOC_OFFSET_SELECT  = "doc_offset_select";
const std::string KEY_DUP  = "dup";
const std::string KEY_DUP_G  = "dup_g";
const std::string KEY_DOCUMENTS  = "documents";
const std::string KEY_P   = "P";
const std::string KEY_P_G   = "P_G";
const std::string KEY_WTP   = "WTP";
const std::string KEY_DOCCNT  = "doccnt";
const std::string KEY_COLLEN  = "collen";
const std::string KEY_DOCBORDER = "docborder";
const std::string KEY_DOCBORDER_RANK = "docborder_rank";
const std::string KEY_DOCBORDER_SELECT = "docborder_select";
const std::string KEY_DOC_LENGTHS = "doclengths";
const std::string KEY_INVFILE_TERM_RANGES = "invfile_term_ranges";
const std::string KEY_INVFILE_PLISTS = "invfile_postings_lists";
const std::string KEY_INVFILE_DOCPERM = "invfile_docperm";
const std::string KEY_INVFILE_IDOCPERM = "invfile_inv_docperm";
const std::string KEY_F_T = "Ft";

const std::string KEY_H = "H";
const std::string KEY_H_SELECT  = "H_select";
const std::string KEY_H_SELECT_0  = "H_select_0";
const std::string KEY_H_SELECT_1  = "H_select_1";
const std::string KEY_H_RANK  = "H_rank";

const std::string KEY_H_LEFT = "H_left";
const std::string KEY_H_LEFT_SELECT_0  = "H_left_select_0";
const std::string KEY_H_LEFT_SELECT_1  = "H_left_select_1";

const std::string KEY_NEXT_OCC = "next_occ";
const std::string KEY_CSA = "csa";
const std::string KEY_MAXTF = "maxtf";
const std::string KEY_MAXDOCLEN = "maxdoclen";
const std::string KEY_QUANTILE_FILTER = "quantile_filter";

const std::string KEY_FILTERED_QUANTILE_FILTER = "filtered_qfilter";
const std::string KEY_FILTERED_QUANTILE_FILTER_RANK = "filtered_qfilter_rank";
const std::string KEY_FILTERED_QUANTILE_FILTER_SELECT = "filtered_qfilter_select";
const std::string KEY_FILTERED_H = "filtered_h";
const std::string KEY_FILTERED_H_SELECT_0 = "filtered_h_select_0";
const std::string KEY_FILTERED_H_SELECT_1 = "filtered_h_select_1";
const std::string KEY_FILTERED_H_RANK = "filtered_h_rank";

const std::string KEY_TAILS = "tails";
const std::string KEY_TAILS_RANK = "tails_rank";
const std::string KEY_TAILS_SELECT = "tails_select";
const std::string KEY_MAXCSTDEPTH = "maxcstdepth";
const std::string KEY_WEIGHTS_RMQ = "weightsrmq";
const std::string KEY_WEIGHTS = "weights";
const std::string KEY_WEIGHTS_G = "weights_g";
const std::string KEY_W_AND_P = "W_and_P";
const std::string KEY_W_AND_P_G = "W_and_P_g";

std::vector<std::string> storage_keys = {
    KEY_DOCCNT,
    KEY_DARRAY,
    KEY_DOCPERM,
    KEY_SADADF,
    KEY_WTD,
    KEY_C,
    KEY_WTC,
    KEY_TMPCST,
    KEY_TMPDUP,
    KEY_DUP,
    KEY_WTDUP,
    KEY_MAXTF,
    KEY_DOCCNT,
    KEY_DOC_LENGTHS,
    KEY_COLLEN,
    KEY_INVFILE_TERM_RANGES,
    KEY_INVFILE_PLISTS,
    KEY_H,
    KEY_U,
    KEY_WTU,
    KEY_CSA,
    KEY_TAILS,
    KEY_NEXT_OCC,
    KEY_WEIGHTS_RMQ,
    sdsl::conf::KEY_TEXT,
    sdsl::conf::KEY_TEXT_INT,
    sdsl::conf::KEY_SA,
    sdsl::conf::KEY_LCP,
    sdsl::conf::KEY_BWT,
    sdsl::conf::KEY_BWT_INT,
    sdsl::conf::KEY_PSI
};

} // end namespace
#endif
