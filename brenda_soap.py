#!/usr/bin/python
#the upper line depends on your operating system
from SOAPpy import WSDL
wsdl = 'http://www.brenda-enzymes.info/soap/brenda.wsdl'
serv = WSDL.Proxy(wsdl)

orgs = ['falciparum','vivax']

for o in orgs:
    results = serv.getKmValue({'organism':'Plasmodium '+o})
    
    
    for r in results:
        ec = r[0][0][1]
        
        print r[0][5][1]
        print r[0][0][1]
        print r[0][4][1]
        print r[0][6][1]
        
        results2 = serv.getRecommendedName({'ecNumber':ec})
        print results2[0][0][1][1]
        break









