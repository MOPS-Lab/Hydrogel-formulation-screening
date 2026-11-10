/* Optional MCMC comparison. Not the production predictor: held-out accuracy
and predictive coverage were worse than the retained GP. */
function createMCMCPredictor(m){
 function cdf(z){const sign=z<0?-1:1,x=Math.abs(z)/Math.sqrt(2),t=1/(1+.3275911*x);const erf=1-(((((1.061405429*t-1.453152027)*t)+1.421413741)*t-.284496736)*t+.254829592)*t*Math.exp(-x*x);return .5*(1+sign*erf);}
 function predict(c,ion,dose){
  const source=ion==='none'?'legacy':'KS',x=[...c,Math.log10(1+dose/.005)].map((v,i)=>(v-m.scaler_mean[i])/m.scaler_scale[i]),n=m.X.length,dist=m.X.map(row=>Math.sqrt(5*row.reduce((s,v,i)=>s+(v-x[i])**2,0))),means=[],variances=[],noise=[];
  for(const d of m.draws){const k=dist.map((v,j)=>{const t=v/d.length;return (1+t+t*t/3)*Math.exp(-t)*(ion===m.ions[j]?1:.3)*(source===m.sources[j]?1:d.rho_source)});let mu=m.y_mean;for(let i=0;i<n;i++)mu+=k[i]*d.alpha[i];const v=new Array(n);for(let i=0;i<n;i++){const offset=i*(i+1)/2;let t=k[i];for(let j=0;j<i;j++)t-=d.L_packed[offset+j]*v[j];v[i]=t/d.L_packed[offset+i];}means.push(mu);variances.push(Math.max(1e-12,1-v.reduce((s,a)=>s+a*a,0)));noise.push(d.noise);}
  const avg=a=>a.reduce((s,v)=>s+v,0)/a.length,mean=avg(means),obsVar=variances.map((v,i)=>v+noise[i]);
  function quant(p){let lo=Math.min(...means.map((v,i)=>v-12*Math.sqrt(obsVar[i]))),hi=Math.max(...means.map((v,i)=>v+12*Math.sqrt(obsVar[i])));for(let it=0;it<65;it++){const mid=(lo+hi)/2,prob=avg(means.map((v,i)=>cdf((mid-v)/Math.sqrt(obsVar[i]))));if(prob<p)lo=mid;else hi=mid;}return 10**((lo+hi)/2);}
  return {G_Pa:10**mean,nominal_predictive80_Pa:[quant(.1),quant(.9)],median_predictive_Pa:quant(.5),variance_log_hyperparameter_means:avg(means.map(v=>(v-mean)**2)),mean_latent_variance_log:avg(variances),mean_noise_variance_log:avg(noise),draws:m.draws.length,production_model:false,conditional_on_gel_formation:true};
 }
 return {predict};
}
