import arch.univariate.base as b, inspect
src = inspect.getsource(b)
for key in ('"robust"', 'def _robust_cov', 'score_cov'):
    i = src.find(key)
    print('===', key, 'at', i)
i = src.find('score_cov')
print(src[max(0, i - 2000):i + 1200])
