# 15418-Final-Project
Parallel A* algorithm using lock-free priority queues.

Summary - Graph search is an important tools in various applications of artificial intelligence from robotic motion planning, last mile optimization to network packet routing. The naive algorithm for graph search being Djikstra's algorithm which is quite computationally inefficient. A* algorithm uses heuristic functions for faster search. We aim to further improve A* algorithm by parallelizing it and achieving a considerable amount of speedup.

Background - A* is an informed, best-first search for finding the minimum cost path on graphs. The algorithm uses a heuristic function which is problem specific which estimates the distance to the goal from a particular vertex to target. The algorithm utilizes an open set of vertices which includes the vertices to search from and a closed set for vertices which have been searched. Starting from a specific starting node of a graph, the algorithm finds the path to the target node having the smallest code by maintaining a tree of paths orignating at the start node and extending those paths one edge every iteration until the target node is reached. At each iteration of the main loop, A* needs to determine which of its paths to extend. It does so based on the cost of the path and an estimate of the cost required to extend the path all the way to the goal. Specifically, A* selects the path that minimizes f(n)=g(n)+h(n) where n is the next node on the path, g(n) is the cost of the path from the start node to n, and h(n) is a heuristic function that estimates the cost of the cheapest path from n to the goal. The heuristic function is problem-specific. If the heuristic function is admissible, meaning that it never overestimates the actual cost to get to the goal, A* is guaranteed to return a least-cost path from start to goal.

A* algorithms use priority queues to perform repeated selections of minimum cost nodes to expand upon which holds the open set. At each step of the algorithm, the node with the lowest f(x) value is removed from the queue, the f and g values of its neighbors are updated accordingly, and these neighbors are added to the queue. When the algorithm reaches the tatgrt node, the f value of that node is then also the cost of the shortest path, since h at the goal is zero in an admissible heuristic. Comdmon heuristic functions include manhattan distance, euclidean distance and diaginal distance.

Initial thoughts on parallelizing A* algorithm include each thread grabbing multiple states from priority queue in order to expand upon multiple paths in parallel. This would lend towards having a global priority queue which is accessed by all threads. In order to maintain correctness, locks for the priority queue data structure would be required. With increased threads, contention for the priority queue could offset the speedup achieved through parallelism. Two approaches to combat this would include using a message passing model in order to distribute work among different threads only allowing the master thread to modify the priority queue. The second more involved mechanism involves implementing a lock-free priority queue. This project is going to implement both methods and compare results alongside the sequential version of the A* algorithm. The project will implemet the lock-free priority queue using a multi-dimensional list based on the paper listed below.

Resources - Bridges Machines:
1. Regular Memory (512x)
○ CPU: 2x AMD EPYC 7742 (2.25-3.40 GHz, 2x64 cores per node)
○ RAM: 256 GB
○ Cache: 256 MB L3 cache, 8 memory channels
○ Local Storage: 3.84 TB NVMe SSD
○ Network: Mellanox ConnectX-6-HDR Infiniband 200Gb/s Adapter

2. V100 GPU (9x)
○ GPU: 8x NVIDIA V100-16GB
○ CPU: 2x Intel Xeon Gold 6148 (2.40-3.70 GHz, 2x20 cores per node)
○ RAM: 192 GB, DDR4-2933
○ Interconnect: PCIe
○ Cache: 2.75 MB last level cache, 6 memory channels
○ Local Storage: 4x 2 TB NVMe SSD
○ Network: 2x Mellanox ConnectX-6-HDR Infiniband 200 Gb/s Adapter

Goals:
- 75% - Implement a parallel A* algorithm using MPI and lock-free priority queue.
- 100% - Improve workload balance of the implementations and perform experiments
- 125% - Implement parallel A* algorithm on GPUs or use a different lock-free priority queue implementation

Schedule-
Week 1- Implement parallel And contention based priority queue
Week 2- Use Message Passing Model to reduce contention
Week 3- Implement lock-free priority queue using multi-dimensional list
Week 4- Optimize workload balance using different techniques
Week 5- Create a distributed set of priority queues for A* implementation on GPUs


Milestone-

1. The main goal of our project is implementing a lock-free priority queue using one the methods which have been published. After much research we were able to find three papers with different implementations of concurrent priority queues which used a mound, skip-list and multidimensional linked list. After reading the three papers carefully we have chosen to work with the multi-dimensional list mainly because it had better experimental performance than the other two data structures. We have not been able to get a working implementation because we were heavily involved in booth before and during carnival. 
2. We have implemented a sequential version of A*, a naive parallel implementation of A* using a locked priority queue. Although we have perfomed conclusive experiments, due to priority queue being mutually exclusive there was not much speedup. Our next goal is to implment a MPI based version where a thread assigns work to nodes to prevent the inefficiency in locking. 

Due to the tremendous booth work and carnival, we were unable to get through to the week 2 deliverables. Because of this our revised schedule looks as follows:
Week 3 - Get a MPI version working and perform first round of experiments
Week 4 - Finish lock free priority queue implementation and perform second round of experiments
Week 5 - Look into work load balance optimizations, parallelization on GPU, approximate concurrent priorirty queue as suggested in comments 

At the demo session, we aim to provide graphs in terms of speedups of different methods, workload balance among nodes and a graph depicting work divided between overhead, synchronization and computation. We also want to provide a visualization for parallel A* graph compared to rudimentary sequential graph search algorithms.
Currently we have no concerns and find ourselves in a good position in terms of understanding the requirements and the project at hand.
