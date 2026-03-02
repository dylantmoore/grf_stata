* test_parity_scope.do -- Validate manifest-based parity scope for gaps 4/5/6
clear all
set more off

local manifest_path "reviews/r_api_manifest.json"
capture confirm file "`manifest_path'"
if _rc {
    local manifest_path "../reviews/r_api_manifest.json"
    capture confirm file "`manifest_path'"
    if _rc {
        display as error "Missing manifest file: reviews/r_api_manifest.json"
        exit 601
    }
}

local in_orthog 0
local in_honesty_enum 0
local in_w_hat 0
local in_y_hat 0
local in_s_hat 0
local in_c_hat 0

local orthog_false 0
local honesty_false 0
local w_hat_true 0
local y_hat_false 0
local s_hat_false 0
local c_hat_false 0

file open fh using "`manifest_path'", read text
file read fh line
while r(eof) == 0 {
    local line_trim = trim(`"`line'"')

    if strpos(`"`line_trim'"', "orthog.boosting") {
        local in_orthog 1
    }
    if strpos(`"`line_trim'"', "honesty.prune.method") {
        local in_honesty_enum 1
    }
    if strpos(`"`line_trim'"', "causal_survival.W.hat") {
        local in_w_hat 1
    }
    if strpos(`"`line_trim'"', "causal_survival.Y.hat") {
        local in_y_hat 1
    }
    if strpos(`"`line_trim'"', "causal_survival.S.hat") {
        local in_s_hat 1
    }
    if strpos(`"`line_trim'"', "causal_survival.C.hat") {
        local in_c_hat 1
    }

    if strpos(`"`line_trim'"', "present") {
        if `in_orthog' & strpos(`"`line_trim'"', "false") {
            local orthog_false 1
            local in_orthog 0
        }
        if `in_honesty_enum' & strpos(`"`line_trim'"', "false") {
            local honesty_false 1
            local in_honesty_enum 0
        }
        if `in_w_hat' & strpos(`"`line_trim'"', "true") {
            local w_hat_true 1
            local in_w_hat 0
        }
        if `in_y_hat' & strpos(`"`line_trim'"', "false") {
            local y_hat_false 1
            local in_y_hat 0
        }
        if `in_s_hat' & strpos(`"`line_trim'"', "false") {
            local s_hat_false 1
            local in_s_hat 0
        }
        if `in_c_hat' & strpos(`"`line_trim'"', "false") {
            local c_hat_false 1
            local in_c_hat 0
        }
    }

    file read fh line
}
file close fh

assert `orthog_false' == 1
assert `honesty_false' == 1
assert `w_hat_true' == 1
assert `y_hat_false' == 1
assert `s_hat_false' == 1
assert `c_hat_false' == 1

display as result "PASS: parity manifest scope checks (gaps 4/5/6)"
