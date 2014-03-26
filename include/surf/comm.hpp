#ifndef SURF_COMM_H
#define SURF_COMM_H

#define REQ_PARSE_ERROR		0
#define REQ_RESPONE_OK		1

#define REQ_TYPE_QRY_OR		0
#define REQ_TYPE_QRY_AND	1
#define REQ_TYPE_QUIT		2

#define REQ_MODE_PROFILE	0
#define REQ_MODE_TIME		1

#define MAX_QRY_LEN		 1024

struct surf_time_resp {
	uint8_t status;
    uint64_t req_id;
    uint64_t qry_id;
    uint64_t qry_len;
    uint64_t k;
    uint64_t result_size;
    uint64_t qry_time;
    uint64_t search_time;
    uint64_t wt_search_space;
    uint64_t postings_evaluated;
    uint64_t postings_total;
};

struct surf_qry_request {
	uint8_t type;
	uint8_t mode;
	uint64_t id;
	uint64_t k;
	char qry_str[MAX_QRY_LEN] = {0};
};


#endif