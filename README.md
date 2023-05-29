# VSS_demo
A practice on verified secret sharing, generally based on Feldman’s verifiable secret sharing protocol.
## For Personal Reference...

The x can be changed to random choice on the finite field, instead of range(k) in my case. x_i just need to be integers falling into the specific range.
----already changed to the version that the x can be chosen by clients agreeing on a integer set in the finite field: see VSS2_main.py.

Also, the initial prime should be of random choice over some large primes, instead my case, which is a fixed number. This prime decides a finite space. And the cost to compute this is one-time overhead.

Two more files are uploaded for full version of VSS and LCC. All operations are done in finite filed of q. Wiki yyds.

Reed Solomon decoding added.

A realization of SNIP is added (reference: Prio: Private, Robust, and Scalable Computation of Aggregate Statistics). SSS and VSS are added in this version of file. And a RS decoding algorithm is also added. 
