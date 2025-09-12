# Versioning

Scientific analysis should always be reproducible. This requires that not only the code is versioned, but also the reference data used (if any). 

BarBeQue is versioned through github, including regular major and minor releases with specific feature sets and software versions.

However, BarBeQue also provides access to several reference databases, some of which are unforunately not versioned. 

Specifically:

| Database | Versioned |
| -------- | --------- |
| All Midori databases | Yes |
| Refseq Mitochondria | Yes |
| Mitofish | No |

Databases like Mitofish are downloaded 'on the day' representing whatever the state of the database is on the day you installed it. This is in contrast to e.g. Midori, which is versioned based on the GenBank release it was built from.

At minimum, we recommend that you make a note of when you have installed the references so you can refer to it properly in any scientific publication you may create later on. You can also make a copy of the reference directoy after you 
have installed all the databases to ensure that you are able to perfectly reproduce it in case your installation gets corrupted down the road. 


