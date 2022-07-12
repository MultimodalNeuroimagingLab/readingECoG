function [ccep_sig,p_th] = ccepPCC_fdr(cross_ccep_p,this_alpha)
%DH 2022
% this_alpha = 0.05;
% function [] = ccep_fdr(cross_ccep_p,this_alpha)
%
[p_sorted] = sort(cross_ccep_p(:));
 
mm = length(p_sorted(~isnan(p_sorted))); % total number of tests, excluding NaNs, those are at the end
 
vector_threshold = this_alpha*[1:length(p_sorted)]/mm;
 
last_sig_val = find((p_sorted'-vector_threshold)<0,1,'last');
 
p_th = vector_threshold(last_sig_val);
ccep_sig = cross_ccep_p<=p_th;