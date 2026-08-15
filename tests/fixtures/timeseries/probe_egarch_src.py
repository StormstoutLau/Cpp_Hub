import arch.univariate.recursions_python as rp
import inspect
src = inspect.getsource(rp.egarch_recursion)
print(src)
