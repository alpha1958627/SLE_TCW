"""Independent ideal AMS exponential oracle, fixed 5000 repetitions per N/D.
Conditional descendants are sampled exactly; no physical lifetime inference.
"""
import json,math,random
def check():
    target=-math.log(.05);report={}
    for D in (1,3):
        rng=random.Random(20261021000+D);estimates=[]
        for repetition in range(5000):
            particles=[rng.expovariate(1.) for i in range(12)];weight=1.
            while True:
                cut=sorted(particles)[D-1]
                if cut>=target:
                    estimates.append(weight*sum(x>=target for x in particles)/12);break
                dead=[i for i,x in enumerate(particles) if x<=cut]
                weight*=1-len(dead)/12
                for i in dead:particles[i]=cut+rng.expovariate(1.)
        mean=math.fsum(estimates)/len(estimates)
        sem=math.sqrt(math.fsum((x-mean)**2 for x in estimates)/(len(estimates)-1)/len(estimates))
        assert abs(mean-.05)<5*sem
        report[str(D)]=dict(N=12,D=D,independent_toy_repetitions=5000,true_P=.05,mean=mean,SEM=sem,
                           nominal_weight_after16=(1-D/12)**16)
    return dict(status='PASS',results=report,scope='ideal continuous AMS conditional-sampling oracle; not q17 physical performance or MFPT')
if __name__=='__main__':print(json.dumps(check(),indent=2))
