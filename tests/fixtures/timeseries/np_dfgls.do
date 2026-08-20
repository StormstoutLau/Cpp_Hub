* np_dfgls.do — Stata dfgls 逐 k MAIC/RMSE 基准生成 (Stata 18 MP, 2026-08-20)
* r(results) 列名 k MAIC SIC RMSE DFGLS; 行序降序 k=kmax..1 (按 k 列值匹配); RMSE=σ̂
version 18.0
set more off

* --- 场景 A: smoke T=200, notrend, maxlag(14) 显式 ---
import delimited using "np_smoke_data.csv", clear case(preserve) stringcols(_all)
destring y, replace force
* dfgls 需 tsset (无时间列 → 用行号; 手册例 lutkepohl2 自带 quarterly tsset)
generate t = _n
tsset t
dfgls y, maxlag(14) notrend
matrix R = r(results)
clear
svmat double R, names(col)
export delimited using "np_stata_baselines.csv", replace

* --- 场景 B: 手册例 lutkepohl2 ln_inv, trend, maxlag(11) 显式 ---
use "https://www.stata-press.com/data/r18/lutkepohl2.dta", clear
* 保留时间变量 qtr (tsset 依赖; dfgls 需 tsset)
keep ln_inv qtr
* lutkepohl2 的 ln_inv 为 float 存储 (~8 位) — recast double (无损) 后全精度
* 导出, 保证 C++ 读入值与 Stata 内部计算值逐位一致 (float→double 精确提升)
recast double ln_inv
format ln_inv %21.0g
export delimited using "lutkepohl2_ln_inv.csv", replace
dfgls ln_inv, maxlag(11)
matrix R2 = r(results)
clear
svmat double R2, names(col)
export delimited using "np_stata_manual_example.csv", replace
exit, clear
