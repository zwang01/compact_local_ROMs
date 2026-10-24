function n = brk_count_fact(reset)
% BRK_COUNT_FACT  Counter of sparse factorizations of the flow saddle-point
% system, the machine-independent cost measure reported in the paper ("fact."
% in the tables; a reduced iteration never factors).
%   brk_count_fact        count one factorization (called before each one)
%   brk_count_fact(0)     reset the counter
%   n = brk_count_fact    count one and return the total
%   n = brk_count_fact([]) return the total without counting
global BRK_NFACT
if isempty(BRK_NFACT), BRK_NFACT = 0; end
if nargin == 0
    BRK_NFACT = BRK_NFACT + 1;
elseif ~isempty(reset)
    BRK_NFACT = reset;
end
n = BRK_NFACT;
end
