# Tests
1. Tried with 2 #pragma omp tasks => It works but it's not that fast (it's actually slower than the sequential version wrt the same input size)
- The problem is that I create an additional task when instead I could just use the thread that creates the first task to directly execute the second "mergesort" call
1. Tried with only 1 #pragma omp tasks => Now it's working better, not as fast as the sequential version, but we have some improvements
1. Tried adding a "final" clause in the "task" construct to stop creating more tasks after a certain number of recursive calls. Now the performance of the parallel version are generally better than the sequential one (well, except for small input sizes, but this is expected given the overhead associated to task creation and scheduling)

# What to write in the pdf
1. Some tables showing the execution times for different input sizes for both the parallel and sequential implementations
    - Maybe some statistics stuff (mean, variance)
1. Performance evaluation for small input sizes
    - Small input => overhead cannot be hidden
    - Maybe we can add an if condition in order to execute the parallel version only when the input size is "large enough"
