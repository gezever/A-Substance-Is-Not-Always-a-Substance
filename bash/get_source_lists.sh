#!/bin/bash

pushd ../data/source
# https://echa.europa.eu/nl/candidate-list-table
curl 'https://echa.europa.eu/nl/candidate-list-table?p_p_id=disslists_WAR_disslistsportlet&p_p_lifecycle=2&p_p_state=normal&p_p_mode=view&p_p_resource_id=exportResults&p_p_cacheability=cacheLevelPage' \
  -X POST \
  -H 'User-Agent: Mozilla/5.0 (X11; Linux x86_64; rv:145.0) Gecko/20100101 Firefox/145.0' \
  -H 'Accept: text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8' \
  -H 'Accept-Language: en-US,en;q=0.5' \
  -H 'Accept-Encoding: gzip, deflate, br, zstd' \
  -H 'Content-Type: application/x-www-form-urlencoded' \
  -H 'Origin: https://echa.europa.eu' \
  -H 'Connection: keep-alive' \
  -H 'Referer: https://echa.europa.eu/nl/candidate-list-table' \
  -H 'Cookie: GUEST_LANGUAGE_ID=nl_NL; COOKIE_SUPPORT=true; cck1=%7B%22cm%22%3Atrue%2C%22all1st%22%3Atrue%7D; _pk_id.c9cba231-694a-48ae-b6c0-eeb7777b02e3.bd53=ddd1ab532f4e552b.1773223190.2.1773230960.1773223190.; LFR_SESSION_STATE_10140=1773233453234; ApplicationGatewayAffinityCORS=6b7407f756da17c61c56c37d83b52516; ApplicationGatewayAffinity=6b7407f756da17c61c56c37d83b52516; JSESSIONID=FEA262BA6709C2A6130A400593C76124.live-2; disclaimer=true' \
  -H 'Upgrade-Insecure-Requests: 1' \
  -H 'Sec-Fetch-Dest: document' \
  -H 'Sec-Fetch-Mode: navigate' \
  -H 'Sec-Fetch-Site: same-origin' \
  -H 'Sec-Fetch-User: ?1' \
  -H 'Priority: u=0, i' \
  --data-raw '_disslists_WAR_disslistsportlet_formDate=1773230956504&_disslists_WAR_disslistsportlet_exportColumns=name%2CecNumber%2CcasNumber%2Chaz_detailed_concern%2Cdte_inclusion%2Cdoc_cat_decision%2Cdoc_cat_iuclid_dossier%2Cdoc_cat_supdoc%2Cdoc_cat_rcom%2Cprc_external_remarks&_disslists_WAR_disslistsportlet_orderByCol=dte_inclusion&_disslists_WAR_disslistsportlet_orderByType=desc&_disslists_WAR_disslistsportlet_searchFormColumns=haz_detailed_concern%2Cdte_inclusion&_disslists_WAR_disslistsportlet_searchFormElements=DROP_DOWN%2CDATE_PICKER&_disslists_WAR_disslistsportlet_substance_identifier_field_key=&_disslists_WAR_disslistsportlet_haz_detailed_concern=&_disslists_WAR_disslistsportlet_dte_inclusionFrom=&_disslists_WAR_disslistsportlet_dte_inclusionTo=&_disslists_WAR_disslistsportlet_total=253&_disslists_WAR_disslistsportlet_exportType=csv'\
  -o candidate-list-of-svhc-for-authorisation-export.csv

# https://chem.echa.europa.eu/obligation-lists/restrictionList
curl 'https://chem.echa.europa.eu/api-obligation-list/v1/restrictionList/fullExport' \
  -H 'User-Agent: Mozilla/5.0 (X11; Linux x86_64; rv:145.0) Gecko/20100101 Firefox/145.0' \
  -H 'Accept: application/json, text/plain, */*' \
  -H 'Accept-Language: en-US,en;q=0.5' \
  -H 'Accept-Encoding: gzip, deflate, br, zstd' \
  -H 'Connection: keep-alive' \
  -H 'Referer: https://chem.echa.europa.eu/obligation-lists/restrictionList' \
  -H 'Cookie: cck1=%7B%22cm%22%3Atrue%2C%22all1st%22%3Atrue%7D; _pk_id.b1c62efe-8fd0-4b63-a615-30d86a21a01e.f839=e4441a60256e9a83.1773223213.2.1773234539.1773223213.; legalNotice=%7B%22accepted%22%3Atrue%2C%22expired%22%3Afalse%2C%22acceptedDate%22%3A1773234539053%7D; _pk_ses.b1c62efe-8fd0-4b63-a615-30d86a21a01e.f839=*' \
  -H 'Sec-Fetch-Dest: empty' \
  -H 'Sec-Fetch-Mode: cors' \
  -H 'Sec-Fetch-Site: same-origin' \
  -H 'Priority: u=0'\
  -o restriction_list_full.xlsx


# https://chem.echa.europa.eu/obligation-lists/candidateList
curl 'https://chem.echa.europa.eu/api-obligation-list/v1/candidateList/fullExport' \
    -H 'User-Agent: Mozilla/5.0 (X11; Linux x86_64; rv:145.0) Gecko/20100101 Firefox/145.0' \
    -H 'Accept: application/json, text/plain, */*' \
    -H 'Accept-Language: en-US,en;q=0.5' \
    -H 'Accept-Encoding: gzip, deflate, br, zstd' \
    -H 'Connection: keep-alive' \
    -H 'Referer: https://chem.echa.europa.eu/obligation-lists/candidateList' \
    -H 'Cookie: cck1=%7B%22cm%22%3Atrue%2C%22all1st%22%3Atrue%7D; _pk_id.b1c62efe-8fd0-4b63-a615-30d86a21a01e.f839=e4441a60256e9a83.1773223213.2.1773234539.1773223213.; legalNotice=%7B%22accepted%22%3Atrue%2C%22expired%22%3Afalse%2C%22acceptedDate%22%3A1773234539053%7D; _pk_ses.b1c62efe-8fd0-4b63-a615-30d86a21a01e.f839=*' \
    -H 'Sec-Fetch-Dest: empty' \
    -H 'Sec-Fetch-Mode: cors' \
    -H 'Sec-Fetch-Site: same-origin' \
    -H 'Priority: u=0'\
    -o candidate_list_full.xlsx

# https://chem.echa.europa.eu/obligation-lists/authorisationList
curl 'https://chem.echa.europa.eu/api-obligation-list/v1/authorisationList/fullExport' \
  -H 'User-Agent: Mozilla/5.0 (X11; Linux x86_64; rv:145.0) Gecko/20100101 Firefox/145.0' \
  -H 'Accept: application/json, text/plain, */*' \
  -H 'Accept-Language: en-US,en;q=0.5' \
  -H 'Accept-Encoding: gzip, deflate, br, zstd' \
  -H 'Connection: keep-alive' \
  -H 'Referer: https://chem.echa.europa.eu/obligation-lists/authorisationList' \
  -H 'Cookie: cck1=%7B%22cm%22%3Atrue%2C%22all1st%22%3Atrue%7D; _pk_id.b1c62efe-8fd0-4b63-a615-30d86a21a01e.f839=e4441a60256e9a83.1773223213.2.1773234539.1773223213.; legalNotice=%7B%22accepted%22%3Atrue%2C%22expired%22%3Afalse%2C%22acceptedDate%22%3A1773234539053%7D; _pk_ses.b1c62efe-8fd0-4b63-a615-30d86a21a01e.f839=*' \
  -H 'Sec-Fetch-Dest: empty' \
  -H 'Sec-Fetch-Mode: cors' \
  -H 'Sec-Fetch-Site: same-origin' \
  -H 'Priority: u=0'\
  -o authorisation_list_full.xlsx

# https://chem.echa.europa.eu/obligation-lists/popsList
curl 'https://chem.echa.europa.eu/api-obligation-list/v1/popsList/fullExport' \
  -H 'User-Agent: Mozilla/5.0 (X11; Linux x86_64; rv:145.0) Gecko/20100101 Firefox/145.0' \
  -H 'Accept: application/json, text/plain, */*' \
  -H 'Accept-Language: en-US,en;q=0.5' \
  -H 'Accept-Encoding: gzip, deflate, br, zstd' \
  -H 'Connection: keep-alive' \
  -H 'Referer: https://chem.echa.europa.eu/obligation-lists/popsList' \
  -H 'Cookie: cck1=%7B%22cm%22%3Atrue%2C%22all1st%22%3Atrue%7D; _pk_id.b1c62efe-8fd0-4b63-a615-30d86a21a01e.f839=e4441a60256e9a83.1773223213.2.1773234539.1773223213.; legalNotice=%7B%22accepted%22%3Atrue%2C%22expired%22%3Afalse%2C%22acceptedDate%22%3A1773234539053%7D; _pk_ses.b1c62efe-8fd0-4b63-a615-30d86a21a01e.f839=*' \
  -H 'Sec-Fetch-Dest: empty' \
  -H 'Sec-Fetch-Mode: cors' \
  -H 'Sec-Fetch-Site: same-origin' \
  -H 'Priority: u=0'\
  -o pops_list_full.xlsx

# https://chem.echa.europa.eu/obligation-lists/euPositiveList
curl 'https://chem.echa.europa.eu/api-obligation-list/v1/euPositiveList/fullExport' \
  -H 'User-Agent: Mozilla/5.0 (X11; Linux x86_64; rv:145.0) Gecko/20100101 Firefox/145.0' \
  -H 'Accept: application/json, text/plain, */*' \
  -H 'Accept-Language: en-US,en;q=0.5' \
  -H 'Accept-Encoding: gzip, deflate, br, zstd' \
  -H 'Connection: keep-alive' \
  -H 'Referer: https://chem.echa.europa.eu/obligation-lists/euPositiveList' \
  -H 'Cookie: cck1=%7B%22cm%22%3Atrue%2C%22all1st%22%3Atrue%7D; _pk_id.b1c62efe-8fd0-4b63-a615-30d86a21a01e.f839=e4441a60256e9a83.1773223213.2.1773234539.1773223213.; legalNotice=%7B%22accepted%22%3Atrue%2C%22expired%22%3Afalse%2C%22acceptedDate%22%3A1773234539053%7D; _pk_ses.b1c62efe-8fd0-4b63-a615-30d86a21a01e.f839=*' \
  -H 'Sec-Fetch-Dest: empty' \
  -H 'Sec-Fetch-Mode: cors' \
  -H 'Sec-Fetch-Site: same-origin' \
  -H 'Priority: u=0'\
  -o eu_positive_list_full.xlsx


# https://chem.echa.europa.eu/obligation-lists/clhList
curl 'https://chem.echa.europa.eu/api-harmonised-list/v1/export?orderBy=indexNumber&orderType=asc&showMembers=false&zoneId=Europe/Amsterdam' \
  -H 'User-Agent: Mozilla/5.0 (X11; Linux x86_64; rv:145.0) Gecko/20100101 Firefox/145.0' \
  -H 'Accept: application/json, text/plain, */*' \
  -H 'Accept-Language: en-US,en;q=0.5' \
  -H 'Accept-Encoding: gzip, deflate, br, zstd' \
  -H 'Connection: keep-alive' \
  -H 'Referer: https://chem.echa.europa.eu/obligation-lists/clhList' \
  -H 'Cookie: cck1=%7B%22cm%22%3Atrue%2C%22all1st%22%3Atrue%7D; _pk_id.b1c62efe-8fd0-4b63-a615-30d86a21a01e.f839=e4441a60256e9a83.1773223213.2.1773234539.1773223213.; legalNotice=%7B%22accepted%22%3Atrue%2C%22expired%22%3Afalse%2C%22acceptedDate%22%3A1773234539053%7D; _pk_ses.b1c62efe-8fd0-4b63-a615-30d86a21a01e.f839=*' \
  -H 'Sec-Fetch-Dest: empty' \
  -H 'Sec-Fetch-Mode: cors' \
  -H 'Sec-Fetch-Site: same-origin' \
  -H 'Priority: u=0'\
  -o Harmonised_List.xlsx


# https://chem.echa.europa.eu/activity-lists/restrictionProcess
curl 'https://chem.echa.europa.eu/api-activity-list/v1/restrictionProcess/fullExport' \
  -H 'User-Agent: Mozilla/5.0 (X11; Linux x86_64; rv:145.0) Gecko/20100101 Firefox/145.0' \
  -H 'Accept: application/json, text/plain, */*' \
  -H 'Accept-Language: en-US,en;q=0.5' \
  -H 'Accept-Encoding: gzip, deflate, br, zstd' \
  -H 'Connection: keep-alive' \
  -H 'Referer: https://chem.echa.europa.eu/activity-lists/restrictionProcess' \
  -H 'Cookie: cck1=%7B%22cm%22%3Atrue%2C%22all1st%22%3Atrue%7D; _pk_id.b1c62efe-8fd0-4b63-a615-30d86a21a01e.f839=e4441a60256e9a83.1773223213.2.1773234539.1773223213.; legalNotice=%7B%22accepted%22%3Atrue%2C%22expired%22%3Afalse%2C%22acceptedDate%22%3A1773234539053%7D; _pk_ses.b1c62efe-8fd0-4b63-a615-30d86a21a01e.f839=*' \
  -H 'Sec-Fetch-Dest: empty' \
  -H 'Sec-Fetch-Mode: cors' \
  -H 'Sec-Fetch-Site: same-origin' \
  -H 'Priority: u=0'\
  -o restriction_process_full.xlsx

# https://chem.echa.europa.eu/activity-lists/svhcIdentification
curl 'https://chem.echa.europa.eu/api-activity-list/v1/svhcIdentification/fullExport' \
  -H 'User-Agent: Mozilla/5.0 (X11; Linux x86_64; rv:145.0) Gecko/20100101 Firefox/145.0' \
  -H 'Accept: application/json, text/plain, */*' \
  -H 'Accept-Language: en-US,en;q=0.5' \
  -H 'Accept-Encoding: gzip, deflate, br, zstd' \
  -H 'Connection: keep-alive' \
  -H 'Referer: https://chem.echa.europa.eu/activity-lists/svhcIdentification' \
  -H 'Cookie: cck1=%7B%22cm%22%3Atrue%2C%22all1st%22%3Atrue%7D; _pk_id.b1c62efe-8fd0-4b63-a615-30d86a21a01e.f839=e4441a60256e9a83.1773223213.2.1773234539.1773223213.; legalNotice=%7B%22accepted%22%3Atrue%2C%22expired%22%3Afalse%2C%22acceptedDate%22%3A1773234539053%7D; _pk_ses.b1c62efe-8fd0-4b63-a615-30d86a21a01e.f839=*' \
  -H 'Sec-Fetch-Dest: empty' \
  -H 'Sec-Fetch-Mode: cors' \
  -H 'Sec-Fetch-Site: same-origin' \
  -H 'Priority: u=0'\
  -o svhc_identification_full.xlsx

# https://chem.echa.europa.eu/activity-lists/authorisationProcess
curl 'https://chem.echa.europa.eu/api-activity-list/v1/authorisationProcess/fullExport' \
  -H 'User-Agent: Mozilla/5.0 (X11; Linux x86_64; rv:145.0) Gecko/20100101 Firefox/145.0' \
  -H 'Accept: application/json, text/plain, */*' \
  -H 'Accept-Language: en-US,en;q=0.5' \
  -H 'Accept-Encoding: gzip, deflate, br, zstd' \
  -H 'Connection: keep-alive' \
  -H 'Referer: https://chem.echa.europa.eu/activity-lists/authorisationProcess' \
  -H 'Cookie: cck1=%7B%22cm%22%3Atrue%2C%22all1st%22%3Atrue%7D; _pk_id.b1c62efe-8fd0-4b63-a615-30d86a21a01e.f839=e4441a60256e9a83.1773223213.2.1773234539.1773223213.; legalNotice=%7B%22accepted%22%3Atrue%2C%22expired%22%3Afalse%2C%22acceptedDate%22%3A1773234539053%7D; _pk_ses.b1c62efe-8fd0-4b63-a615-30d86a21a01e.f839=*' \
  -H 'Sec-Fetch-Dest: empty' \
  -H 'Sec-Fetch-Mode: cors' \
  -H 'Sec-Fetch-Site: same-origin' \
  -H 'Priority: u=0'\
  -o authorisation_process_full.xlsx


# https://chem.echa.europa.eu/activity-lists/dossierEvaluation
curl 'https://chem.echa.europa.eu/api-activity-list/v1/dossierEvaluation/fullExport' \
  -H 'User-Agent: Mozilla/5.0 (X11; Linux x86_64; rv:145.0) Gecko/20100101 Firefox/145.0' \
  -H 'Accept: application/json, text/plain, */*' \
  -H 'Accept-Language: en-US,en;q=0.5' \
  -H 'Accept-Encoding: gzip, deflate, br, zstd' \
  -H 'Connection: keep-alive' \
  -H 'Referer: https://chem.echa.europa.eu/activity-lists/dossierEvaluation' \
  -H 'Cookie: cck1=%7B%22cm%22%3Atrue%2C%22all1st%22%3Atrue%7D; _pk_id.b1c62efe-8fd0-4b63-a615-30d86a21a01e.f839=e4441a60256e9a83.1773223213.2.1773234539.1773223213.; legalNotice=%7B%22accepted%22%3Atrue%2C%22expired%22%3Afalse%2C%22acceptedDate%22%3A1773234539053%7D; _pk_ses.b1c62efe-8fd0-4b63-a615-30d86a21a01e.f839=*' \
  -H 'Sec-Fetch-Dest: empty' \
  -H 'Sec-Fetch-Mode: cors' \
  -H 'Sec-Fetch-Site: same-origin' \
  -H 'Priority: u=0'\
  -o dossier_evaluation_full.xlsx


# https://chem.echa.europa.eu/activity-lists/clhProcess
curl 'https://chem.echa.europa.eu/api-activity-list/v1/clhProcess/fullExport' \
  -H 'User-Agent: Mozilla/5.0 (X11; Linux x86_64; rv:145.0) Gecko/20100101 Firefox/145.0' \
  -H 'Accept: application/json, text/plain, */*' \
  -H 'Accept-Language: en-US,en;q=0.5' \
  -H 'Accept-Encoding: gzip, deflate, br, zstd' \
  -H 'Connection: keep-alive' \
  -H 'Referer: https://chem.echa.europa.eu/activity-lists/clhProcess' \
  -H 'Cookie: cck1=%7B%22cm%22%3Atrue%2C%22all1st%22%3Atrue%7D; _pk_id.b1c62efe-8fd0-4b63-a615-30d86a21a01e.f839=e4441a60256e9a83.1773223213.2.1773234539.1773223213.; legalNotice=%7B%22accepted%22%3Atrue%2C%22expired%22%3Afalse%2C%22acceptedDate%22%3A1773234539053%7D; _pk_ses.b1c62efe-8fd0-4b63-a615-30d86a21a01e.f839=*' \
  -H 'Sec-Fetch-Dest: empty' \
  -H 'Sec-Fetch-Mode: cors' \
  -H 'Sec-Fetch-Site: same-origin' \
  -H 'Priority: u=0'\
  -o clh_process_full.xlsx


# https://chem.echa.europa.eu/activity-lists/substanceEvaluation
curl 'https://chem.echa.europa.eu/api-activity-list/v1/substanceEvaluation/fullExport' \
  -H 'User-Agent: Mozilla/5.0 (X11; Linux x86_64; rv:145.0) Gecko/20100101 Firefox/145.0' \
  -H 'Accept: application/json, text/plain, */*' \
  -H 'Accept-Language: en-US,en;q=0.5' \
  -H 'Accept-Encoding: gzip, deflate, br, zstd' \
  -H 'Connection: keep-alive' \
  -H 'Referer: https://chem.echa.europa.eu/activity-lists/substanceEvaluation' \
  -H 'Cookie: cck1=%7B%22cm%22%3Atrue%2C%22all1st%22%3Atrue%7D; _pk_id.b1c62efe-8fd0-4b63-a615-30d86a21a01e.f839=e4441a60256e9a83.1773223213.2.1773234539.1773223213.; legalNotice=%7B%22accepted%22%3Atrue%2C%22expired%22%3Afalse%2C%22acceptedDate%22%3A1773234539053%7D; _pk_ses.b1c62efe-8fd0-4b63-a615-30d86a21a01e.f839=*' \
  -H 'Sec-Fetch-Dest: empty' \
  -H 'Sec-Fetch-Mode: cors' \
  -H 'Sec-Fetch-Site: same-origin' \
  -H 'Priority: u=0'\
  -o substance_evaluation_full.xlsx

# https://chem.echa.europa.eu/activity-lists/popsProcess
curl 'https://chem.echa.europa.eu/api-activity-list/v1/popsProcess/fullExport' \
  -H 'User-Agent: Mozilla/5.0 (X11; Linux x86_64; rv:145.0) Gecko/20100101 Firefox/145.0' \
  -H 'Accept: application/json, text/plain, */*' \
  -H 'Accept-Language: en-US,en;q=0.5' \
  -H 'Accept-Encoding: gzip, deflate, br, zstd' \
  -H 'Connection: keep-alive' \
  -H 'Referer: https://chem.echa.europa.eu/activity-lists/popsProcess' \
  -H 'Cookie: cck1=%7B%22cm%22%3Atrue%2C%22all1st%22%3Atrue%7D; _pk_id.b1c62efe-8fd0-4b63-a615-30d86a21a01e.f839=e4441a60256e9a83.1773223213.2.1773234539.1773223213.; legalNotice=%7B%22accepted%22%3Atrue%2C%22expired%22%3Afalse%2C%22acceptedDate%22%3A1773234539053%7D; _pk_ses.b1c62efe-8fd0-4b63-a615-30d86a21a01e.f839=*' \
  -H 'Sec-Fetch-Dest: empty' \
  -H 'Sec-Fetch-Mode: cors' \
  -H 'Sec-Fetch-Site: same-origin' \
  -H 'Priority: u=0'\
  -o pops_process_full.xlsx

#   https://chem.echa.europa.eu/activity-lists/pbtAssessment
curl 'https://chem.echa.europa.eu/api-activity-list/v1/pbtAssessment/export?orderBy=currentStageDate&orderType=desc&showMembers=false&zoneId=Europe/Amsterdam' \
  -H 'User-Agent: Mozilla/5.0 (X11; Linux x86_64; rv:145.0) Gecko/20100101 Firefox/145.0' \
  -H 'Accept: application/json, text/plain, */*' \
  -H 'Accept-Language: en-US,en;q=0.5' \
  -H 'Accept-Encoding: gzip, deflate, br, zstd' \
  -H 'Connection: keep-alive' \
  -H 'Referer: https://chem.echa.europa.eu/activity-lists/pbtAssessment' \
  -H 'Cookie: cck1=%7B%22cm%22%3Atrue%2C%22all1st%22%3Atrue%7D; _pk_id.b1c62efe-8fd0-4b63-a615-30d86a21a01e.f839=e4441a60256e9a83.1773223213.2.1773234539.1773223213.; legalNotice=%7B%22accepted%22%3Atrue%2C%22expired%22%3Afalse%2C%22acceptedDate%22%3A1773234539053%7D; _pk_ses.b1c62efe-8fd0-4b63-a615-30d86a21a01e.f839=*' \
  -H 'Sec-Fetch-Dest: empty' \
  -H 'Sec-Fetch-Mode: cors' \
  -H 'Sec-Fetch-Site: same-origin' \
  -H 'Priority: u=0'\
  -o pbt_assessment.xlsx

#   https://chem.echa.europa.eu/activity-lists/edAssessment
curl 'https://chem.echa.europa.eu/api-activity-list/v1/edAssessment/export?orderBy=currentStageDate&orderType=desc&showMembers=false&zoneId=Europe/Amsterdam' \
  -H 'User-Agent: Mozilla/5.0 (X11; Linux x86_64; rv:145.0) Gecko/20100101 Firefox/145.0' \
  -H 'Accept: application/json, text/plain, */*' \
  -H 'Accept-Language: en-US,en;q=0.5' \
  -H 'Accept-Encoding: gzip, deflate, br, zstd' \
  -H 'Connection: keep-alive' \
  -H 'Referer: https://chem.echa.europa.eu/activity-lists/edAssessment' \
  -H 'Cookie: cck1=%7B%22cm%22%3Atrue%2C%22all1st%22%3Atrue%7D; _pk_id.b1c62efe-8fd0-4b63-a615-30d86a21a01e.f839=e4441a60256e9a83.1773223213.2.1773234539.1773223213.; legalNotice=%7B%22accepted%22%3Atrue%2C%22expired%22%3Afalse%2C%22acceptedDate%22%3A1773234539053%7D; _pk_ses.b1c62efe-8fd0-4b63-a615-30d86a21a01e.f839=*' \
  -H 'Sec-Fetch-Dest: empty' \
  -H 'Sec-Fetch-Mode: cors' \
  -H 'Sec-Fetch-Site: same-origin' \
  -H 'Priority: u=0'\
  -o ed_assessment.xlsx


# Show all substances https://chem.echa.europa.eu/
curl https://chem.echa.europa.eu/api-substance/v1/substance/generated-export -o reach_registrations.xlsx

# OSPAR List of Chemicals for Priority Action
curl 'https://www.ospar.org/documents?d=32745' \
  -H 'User-Agent: Mozilla/5.0 (X11; Linux x86_64; rv:145.0) Gecko/20100101 Firefox/145.0' \
  -H 'Accept: text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8' \
  -H 'Accept-Language: en-US,en;q=0.5' \
  -H 'Accept-Encoding: gzip, deflate, br, zstd' \
  -H 'Referer: https://www.ospar.org/work-areas/hasec/hazardous-substances/priority-action' \
  -H 'Connection: keep-alive' \
  -H 'Cookie: cf_clearance=YImSPclZAyUwwJ.WnLot8zZvpJtfE1AhSMvXhv7oDMY-1773237490-1.2.1.1-xBZV5YWfg6dB0giqAmbQ1rC0pm62d_muTAkdsuGGT_83TmZw9mHUPO8P0DmQf34hq8zCMMBzQ1g6Fe.SPBdPTFUL7Eq7AmSZosb65p7beJxAqfq1G9DEVgTHNgRXs5QBDr.U3DHSrSb0uCD9SaDo0AUrUw1a5WMt9pItvj06UlR83mQGfkLoj.sXkLt._OwKCkhR3wvkkVFtHAgjzsMGjwG_vm32SZPR8O6L91LMNnU; _ga_ZS9G6HV5C8=GS2.1.s1773237280$o1$g1$t1773237569$j52$l0$h0; _ga=GA1.1.654486601.1773237281; _ga_436QV5PZQJ=GS2.1.s1773237280$o1$g1$t1773237561$j60$l0$h0; _gid=GA1.2.2075219325.1773237281; wires=598ad0cf7ae6a3d610c3f1e36f8a01eb' \
  -H 'Upgrade-Insecure-Requests: 1' \
  -H 'Sec-Fetch-Dest: document' \
  -H 'Sec-Fetch-Mode: navigate' \
  -H 'Sec-Fetch-Site: same-origin' \
  -H 'Sec-Fetch-User: ?1' \
  -H 'Priority: u=0, i' \
  -H 'TE: trailers'\
  -o 04-12e_agreement_list_of_chemicals_for_priority_action_update_2023-24.pdf

# https://www.efsa.europa.eu/en/applications/pesticides
curl 'https://ec.europa.eu/food/plant/pesticides/eu-pesticides-database/backend/api/active_substance/export' \
  --compressed \
  -X POST \
  -H 'User-Agent: Mozilla/5.0 (X11; Linux x86_64; rv:145.0) Gecko/20100101 Firefox/145.0' \
  -H 'Accept: application/json, text/plain, */*' \
  -H 'Accept-Language: en-US,en;q=0.5' \
  -H 'Accept-Encoding: gzip, deflate, br, zstd' \
  -H 'X-Requested-With: XMLHttpRequest' \
  -H 'Cache-Control: No-Cache' \
  -H 'Content-Type: application/json' \
  -H 'Origin: https://ec.europa.eu' \
  -H 'Connection: keep-alive' \
  -H 'Referer: https://ec.europa.eu/food/plant/pesticides/eu-pesticides-database/start/screen/active-substances' \
  -H 'Cookie: _pk_id.e1483677-4df3-4239-a9ad-9c01253e73e0.5033=60d8f6e545723136.1773238004.1.1773238004..; _pk_ses.e1483677-4df3-4239-a9ad-9c01253e73e0.5033=*; cck1=%7B%22cm%22%3Atrue%2C%22all1st%22%3Atrue%7D' \
  -H 'Sec-Fetch-Dest: empty' \
  -H 'Sec-Fetch-Mode: cors' \
  -H 'Sec-Fetch-Site: same-origin' \
  -H 'Priority: u=0' \
  --data-raw '["1062","351","1329","357","358","359","360","1330","361","696","352","362","355","349","1499","1241","350","1328","363","364","1331","1340","1240","365","366","367","1017","368","810","1341","369","370","371","372","1342","1496","373","374","1245","1246","375","376","1512","377","1343","1380","353","1226","356","1500","354","378","379","382","825","1344","383","381","1432","384","385","386","1203","387","826","388","389","391","380","1561","390","827","393","394","395","617","396","397","392","828","398","399","400","401","402","403","404","405","406","407","408","1032","1013","1033","1337","855","1034","1035","1037","1038","1039","1040","1042","1041","1043","991","1044","1045","1046","1259","1307","1505","1048","1049","1209","1050","1051","1052","1418","1053","1014","1054","1055","1056","1239","863","1299","309","310","265","311","312","313","314","315","316","317","319","320","321","322","1424","1457","323","324","325","326","1015","1280","327","328","705","329","1242","330","331","334","865","335","336","337","1569","339","340","341","342","343","345","346","347","348","1016","409","1217","1506","1573","1260","410","1236","872","411","413","1401","414","1549","415","416","417","418","419","420","421","422","423","424","425","426","427","428","1558","1018","1257","1476","1333","1197","1448","1198","1550","1078","992","1262","1446","1079","429","1261","1278","430","1264","1502","1269","1301","861","1270","1271","1272","1273","1463","433","1571","1265","1570","1483","1555","1514","434","1510","435","436","437","438","1183","1336","1275","1305","1494","1339","1282","1184","1281","1511","1215","439","1415","1524","440","1345","441","442","443","444","445","446","1346","447","1019","449","450","1347","451","452","453","454","455","1204","1348","456","457","458","459","1266","1243","1349","866","460","461","462","463","1350","464","465","466","1393","467","675","1420","1495","1436","1391","468","469","470","471","472","473","474","475","476","477","478","479","480","481","482","483","484","1381","485","486","487","682","488","489","490","829","1185","491","492","493","494","495","830","496","497","498","499","500","1438","831","501","502","503","1277","1508","504","505","832","506","507","508","509","1454","510","511","1456","512","513","514","515","1068","516","517","518","1065","519","520","521","1490","1193","1072","522","523","525","524","1351","526","527","528","529","530","1352","531","532","987","533","534","1353","1020","535","536","537","538","1220","662","1354","539","540","704","541","542","1355","543","544","545","546","547","1356","548","549","550","551","552","1357","553","554","1477","1021","1213","1279","555","1421","1551","556","558","559","560","1225","561","562","563","564","565","981","766","566","867","1449","568","569","570","1221","571","580","572","573","1022","574","852","575","1255","576","577","1080","1482","578","579","581","582","1382","1083","583","584","1289","585","586","587","588","679","833","978","589","590","591","592","593","834","594","1023","595","596","1515","597","600","693","598","599","835","602","603","604","836","605","606","607","610","608","609","611","612","613","614","1358","615","616","618","1359","619","620","621","622","623","1360","624","625","626","1190","627","1361","628","629","630","631","632","1362","633","634","635","1292","636","637","1363","638","639","640","641","642","643","1214","644","1428","988","645","646","868","1491","647","648","649","650","651","652","653","1394","654","655","656","657","658","659","1383","690","1025","1365","691","680","692","1064","837","694","695","1366","1451","698","1368","702","697","1442","1422","838","699","1367","700","701","106","84","677","85","1077","676","812","1210","1311","86","87","88","89","90","999","813","91","92","93","94","95","96","814","1066","97","663","98","99","100","101","102","103","104","301","107","108","1191","118","109","110","1000","111","112","846","114","116","117","119","120","121","122","123","124","125","1194","2","815","4","5","6","7","8","816","9","10","11","1286","12","13","817","683","14","1474","15","16","1300","17","18","19","818","669","20","21","22","23","1310","1002","24","25","819","26","27","28","1414","29","1472","30","31","32","1312","33","34","660","35","36","37","1313","38","1371","1263","717","718","719","670","720","721","722","723","724","725","1003","1444","726","727","1186","728","729","731","732","733","734","735","736","737","738","1059","739","740","989","741","742","743","744","745","746","747","748","1004","749","750","751","752","753","1400","754","755","756","757","758","759","1283","1507","1211","760","761","762","763","1036","764","765","79","1314","811","1252","80","81","769","1413","82","83","767","768","771","1005","772","1372","773","770","774","126","1498","127","128","1492","129","130","1001","131","132","133","134","135","1572","136","1315","137","138","986","139","140","1316","1364","687","851","703","688","689","993","1061","227","228","229","1267","230","231","232","233","1293","979","234","235","236","237","938","238","239","240","241","242","1199","1233","243","244","245","1189","246","1024","671","665","247","248","249","250","251","1488","252","253","254","255","256","257","412","1481","1291","995","258","260","1303","261","262","263","264","1223","266","1208","1519","1441","267","666","268","269","1379","1323","270","271","1452","1324","279","891","1467","272","273","820","274","275","276","277","1200","996","278","821","1244","280","281","282","822","1332","283","284","285","286","287","823","1326","288","289","203","204","205","1318","206","207","208","980","209","1479","1256","1319","1464","1465","1466","210","211","212","213","1320","214","215","216","217","218","843","219","1009","1321","1523","220","221","222","844","847","223","226","224","1322","225","1385","1212","876","1390","1069","877","878","1195","845","1306","879","880","1232","881","882","1287","1288","883","884","848","885","886","887","888","889","890","892","1416","893","896","914","895","897","898","899","900","1258","901","1513","902","903","904","905","1386","906","907","908","909","1387","910","911","912","913","915","916","917","918","919","920","921","922","1304","1458","1453","39","923","1425","924","1070","925","926","927","928","929","930","931","1290","932","1076","933","934","935","1475","936","1431","937","939","940","1517","1417","941","942","1010","943","944","945","946","947","948","949","950","951","952","953","954","1309","958","1075","959","960","1192","961","853","962","963","964","1230","982","983","1334","1335","1187","1011","965","966","967","968","869","969","970","971","972","973","1388","974","849","975","1294","1295","1296","1426","1427","976","977","1389","850","1029","854","1030","1031","1074","290","291","824","292","672","684","997","293","294","295","296","298","299","300","302","303","304","305","306","307","308","40","41","42","43","44","46","48","49","1480","50","51","1201","1302","52","53","54","55","56","661","57","58","1073","59","60","61","62","63","664","64","65","998","66","67","68","69","70","71","72","73","74","75","76","77","78","706","707","708","1369","984","709","710","711","712","1370","713","714","715","716","1084","1085","1285","864","1086","1308","1405","1087","1219","1088","1089","1090","1406","1091","1092","1093","1094","1095","1096","1407","1097","1012","681","1098","1099","985","862","1100","1248","1408","1067","1101","1468","1102","1104","1533","1103","1106","1409","1107","1108","1047","1109","1110","1410","1111","1112","1113","1114","870","1115","856","1120","1116","1117","1026","1118","1119","1188","1121","857","1063","873","1122","1123","1124","1137","1202","1196","667","874","1276","1254","1125","1126","858","1127","1128","1392","1129","1130","1131","1132","859","1133","1134","1136","668","1027","860","1138","1139","1140","1141","1142","1143","1144","1145","1146","1147","1148","1235","1149","1150","1504","1151","1152","1153","1154","1155","1157","1159","1395","1156","1158","1135","1160","1161","1162","1163","1164","1165","1166","1167","1168","1169","1170","1171","1423","1172","1173","1530","1528","1529","1552","1081","1411","1175","1176","1206","1177","673","1412","1179","1058","1180","1028","1181","1182","45","1234","777","778","685","1419","875","775","776","779","780","1374","781","782","839","783","784","1218","785","1375","786","840","787","788","789","1060","1376","790","1384","791","792","793","1377","794","795","841","796","797","798","799","800","842","1378","801","802","803","804","805","806","807","808","809","141","686","1082","1071","142","143","144","1284","145","146","147","148","149","150","1535","154","152","153","155","1006","156","157","158","159","1222","160","161","162","163","164","1403","169","1484","1396","165","1397","1398","674","1402","1297","1298","1433","1231","1268","167","1205","168","1455","1216","170","171","172","173","174","175","176","1249","1478","177","1007","178","179","180","181","182","183","184","185","186","1575","994","187","1224","188","189","190","191","192","1338","193","1207","1440","194","195","196","1399","1429","197","1439","199","1317","1008","200","201","202","557","601","1473","259","1253","955","956","957","1404","1237","1238","1373","151","198"]'\
  -o Pesticides_ActiveSubstanceExport.xlsx


popd
