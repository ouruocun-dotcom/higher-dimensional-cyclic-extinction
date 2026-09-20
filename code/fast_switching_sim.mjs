import fs from 'node:fs';
import path from 'node:path';

function arg(name, fallback) {
  const i = process.argv.indexOf(`--${name}`);
  return i >= 0 ? process.argv[i + 1] : fallback;
}
const cfg = {
  q: +arg('q', 3), trials: +arg('trials', 1000), kMean: +arg('k-mean', 280),
  deltaK: +arg('delta-k', 0.7), s: +arg('s', 0.5), nu: +arg('nu', 20),
  heterogeneity: +arg('heterogeneity', 0.55), phase: +arg('phase', 0),
  seed: +arg('seed', 20260919), output: arg('output', 'research/fast_switching_node.json'),
  fullFixation: arg('full-fixation', 'false') === 'true', topology: arg('topology', 'ring'),
  design: arg('design', 'all'), leakage: +arg('leakage', 0)
};

function mulberry32(seed) {
  return function () {
    let t = seed += 0x6D2B79F5;
    t = Math.imul(t ^ t >>> 15, t | 1);
    t ^= t + Math.imul(t ^ t >>> 7, t | 61);
    return ((t ^ t >>> 14) >>> 0) / 4294967296;
  };
}

function matrix(q, h, phase, topology, leakage=0) {
  const w = Array.from({length:q}, (_,i) => Math.exp(h*Math.cos(2*Math.PI*i/q+phase)));
  const mean = w.reduce((a,b)=>a+b,0)/q;
  const A = Array.from({length:q},()=>new Float64Array(q));
  if(topology==='ring'){
    for (let i=0;i<q;i++) { const j=(i+1)%q, z=w[i]/mean; A[i][j]=z; A[j][i]=-z; }
  } else if(topology==='tournament'){
    if(q%2===0)throw Error('cyclic tournament requires odd q');
    const m=(q-1)/2;
    for(let i=0;i<q;i++)for(let d=1;d<=m;d++){
      const j=(i+d)%q,z=Math.sqrt(w[i]*w[j])/mean/m;A[i][j]=z;A[j][i]=-z;
    }
  } else if(topology==='leaky_ring'){
    if(q%2===0)throw Error('leaky ring requires odd q');
    for (let i=0;i<q;i++) { const j=(i+1)%q, z=w[i]/mean; A[i][j]=z; A[j][i]=-z; }
    const m=(q-1)/2;
    for(let i=0;i<q;i++)for(let d=1;d<=m;d++){
      const j=(i+d)%q;
      if(Math.abs(A[i][j])<1e-15){const z=leakage*Math.sqrt(w[i]*w[j])/mean/m;A[i][j]=z;A[j][i]=-z;}
    }
  } else throw Error(`unknown topology ${topology}`);
  return A;
}

function nullVector(A) {
  const q=A.length, n=q-1, m=Array.from({length:n},(_,i)=>{
    const r=new Float64Array(n+1); for(let j=0;j<n;j++)r[j]=A[i][j]; r[n]=-A[i][q-1]; return r;
  });
  for(let c=0;c<n;c++){
    let p=c; for(let r=c+1;r<n;r++) if(Math.abs(m[r][c])>Math.abs(m[p][c]))p=r;
    [m[c],m[p]]=[m[p],m[c]]; const d=m[c][c]; if(Math.abs(d)<1e-12)throw Error('singular');
    for(let j=c;j<=n;j++)m[c][j]/=d;
    for(let r=0;r<n;r++)if(r!==c){const f=m[r][c];for(let j=c;j<=n;j++)m[r][j]-=f*m[c][j];}
  }
  const x=Array(q).fill(1); for(let i=0;i<n;i++)x[i]=m[i][n];
  if(x.reduce((a,b)=>a+b,0)<0)for(let i=0;i<q;i++)x[i]*=-1;
  if(x.some(v=>v<=0))throw Error('nonpositive equilibrium');
  const z=x.reduce((a,b)=>a+b,0); return x.map(v=>v/z);
}

function initialCounts(total,x){
  const raw=x.map(v=>v*total), n=raw.map(v=>Math.max(1,Math.floor(v)));
  let used=n.reduce((a,b)=>a+b,0);
  while(used<total){let b=0;for(let i=1;i<n.length;i++)if(raw[i]-n[i]>raw[b]-n[b])b=i;n[b]++;used++;}
  while(used>total){let b=-1;for(let i=0;i<n.length;i++)if(n[i]>1&&(b<0||n[i]-raw[i]>n[b]-raw[b]))b=i;n[b]--;used--;}
  return n;
}

function trajectory(rng,A,xstar,s,Kmean,delta,nu,switching,initialTotal,fullFixation){
  const q=A.length,n=initialCounts(initialTotal,xstar),birth=new Float64Array(q),death=new Float64Array(q);
  let total=initialTotal,env=rng()<.5?-1:1,time=0,events=0,alive=q;
  const order=[],times=[];
  const targetAlive=fullFixation?1:q-1;
  while(alive>targetAlive){
    events++; let bs=0,ds=0;
    for(let i=0;i<q;i++){
      if(n[i]===0){birth[i]=death[i]=0;continue;}
      let payoff=0;for(let j=0;j<q;j++)payoff+=A[i][j]*n[j]/total;
      const fit=1+s*payoff;if(fit<=0)throw Error('nonpositive fitness'); birth[i]=n[i]*fit;bs+=birth[i];
    }
    const K=switching?Kmean*(1+delta*env):Kmean;
    for(let i=0;i<q;i++){death[i]=n[i]*total/K;ds+=death[i];}
    const sr=switching?nu:0,rate=bs+ds+sr; time+=-Math.log(Math.max(rng(),Number.MIN_VALUE))/rate;
    const u=rng()*rate;
    if(u<bs){let c=0;for(let i=0;i<q;i++)if((c+=birth[i])>=u){n[i]++;total++;break;}}
    else if(u<bs+ds){let c=0,t=u-bs;for(let i=0;i<q;i++)if((c+=death[i])>=t){n[i]--;total--;if(n[i]===0){order.push(i);times.push(time);alive--;}break;}}
    else env=-env;
  }
  return {extinct:order[0],time:times[0],events,order,times};
}

const A=matrix(cfg.q,cfg.heterogeneity,cfg.phase,cfg.topology,cfg.leakage),xstar=nullVector(A);
const kEff=Math.max(cfg.q,Math.round(cfg.kMean*(1-cfg.deltaK**2))),sEff=cfg.s*(1-cfg.deltaK**2);
const allDesigns=[
  ['switching',cfg.s,cfg.kMean,true,kEff],['static_naive',cfg.s,cfg.kMean,false,cfg.kMean],
  ['static_effective_s',sEff,cfg.kMean,false,cfg.kMean],['static_effective_k',cfg.s,kEff,false,kEff]
];
const designs=cfg.design==='all'?allDesigns:allDesigns.filter(d=>d[0]===cfg.design);
if(designs.length===0)throw Error(`unknown design ${cfg.design}`);
const results={};
for(let d=0;d<designs.length;d++){
  const [name,s,K,sw,init]=designs[d],counts=Array(cfg.q).fill(0);
  const stageCounts=Array.from({length:cfg.q-1},()=>Array(cfg.q).fill(0));
  const stageAbsoluteTime=Array(cfg.q-1).fill(0),stageResidenceTime=Array(cfg.q-1).fill(0);
  const stageResidenceTimeSq=Array(cfg.q-1).fill(0);
  const orders=new Map(),finalPairs=new Map();let tsum=0,esum=0;
  for(let t=0;t<cfg.trials;t++){
    const rng=mulberry32((cfg.seed+Math.imul(1000003,d+1)+Math.imul(9176,t+1))>>>0);
    const r=trajectory(rng,A,xstar,s,K,cfg.deltaK,cfg.nu,sw,init,cfg.fullFixation);counts[r.extinct]++;tsum+=r.time;esum+=r.events;
    if(cfg.fullFixation){
      for(let stage=0;stage<r.order.length;stage++){
        stageCounts[stage][r.order[stage]]++;
        stageAbsoluteTime[stage]+=r.times[stage];
        const residence=r.times[stage]-(stage===0?0:r.times[stage-1]);
        stageResidenceTime[stage]+=residence;
        stageResidenceTimeSq[stage]+=residence*residence;
      }
      const key=r.order.join('-');orders.set(key,(orders.get(key)||0)+1);
      const winner=Array.from({length:cfg.q},(_,i)=>i).find(i=>!r.order.includes(i));
      const lastExtinct=r.order[r.order.length-1];
      const pair=[winner,lastExtinct].sort((a,b)=>a-b),pairKey=pair.join('-');
      const residence=r.times[r.times.length-1]-(r.times.length===1?0:r.times[r.times.length-2]);
      const pairRow=finalPairs.get(pairKey)||{count:0,sumResidenceTime:0,sumSqResidenceTime:0,interactionMagnitude:Math.abs(A[pair[0]][pair[1]])};
      pairRow.count++;pairRow.sumResidenceTime+=residence;pairRow.sumSqResidenceTime+=residence*residence;
      finalPairs.set(pairKey,pairRow);
    }
  }
  const orderCounts=Object.fromEntries([...orders.entries()].sort((a,b)=>b[1]-a[1]));
  const finalPairStats=Object.fromEntries([...finalPairs.entries()].sort((a,b)=>b[1].count-a[1].count).map(([key,row])=>[key,{
    ...row,meanResidenceTime:row.sumResidenceTime/row.count
  }]));
  results[name]={counts,probabilities:counts.map(v=>v/cfg.trials),meanFirstExtinctionTime:tsum/cfg.trials,meanEvents:esum/cfg.trials,
    stageCounts,stageProbabilities:stageCounts.map(row=>row.map(v=>v/cfg.trials)),
    meanStageAbsoluteTime:stageAbsoluteTime.map(v=>v/cfg.trials),
    meanStageResidenceTime:stageResidenceTime.map(v=>v/cfg.trials),
    seStageResidenceTime:stageResidenceTime.map((v,i)=>{
      const mean=v/cfg.trials;
      const variance=Math.max(0,(stageResidenceTimeSq[i]-cfg.trials*mean*mean)/(cfg.trials-1));
      return Math.sqrt(variance/cfg.trials);
    }),orderCounts,finalPairStats};
}
const out={config:cfg,derived:{kEff,sEff},equilibrium:xstar,designs:results};
fs.mkdirSync(path.dirname(cfg.output),{recursive:true});fs.writeFileSync(cfg.output,JSON.stringify(out,null,2));
console.log(JSON.stringify({output:cfg.output,derived:out.derived,meanEvents:Object.fromEntries(Object.entries(results).map(([k,v])=>[k,v.meanEvents]))}));
