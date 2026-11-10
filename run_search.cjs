const fs=require('fs'),vm=require('vm'),path=require('path');
const R=path.join(__dirname,'analysis_source'),out=path.join(__dirname,'revision_analysis');fs.mkdirSync(out,{recursive:true});
const load=x=>JSON.parse(fs.readFileSync(path.join(R,x),'utf8'));
const ctx={console,setTimeout};vm.createContext(ctx);
vm.runInContext(fs.readFileSync(R+'/model/hydrogel.js','utf8')+'\n'+fs.readFileSync(R+'/model/monte_carlo.js','utf8'),ctx);
const meta=load('model/metadata.json'),model=ctx.createHydrogelModel(load('model/joint_gp.json'),meta,load('model/diagnostics.json'),load('model/alternative.json')),mc=ctx.createMonteCarlo(model,meta);
(async()=>{
 const results=[];
 for(const target of [5000,8000,10000,12000,15000,20000,30000]){
  const r=await mc.robustInverse({target_Pa:target,tolerance_percent:10,family:'ALG+TOCNF',ion:'Ca',supported_only:true},{n:5000,seed:20261014});
  results.push(r);fs.writeFileSync(out+'/search_results.json',JSON.stringify(results,null,2));
  console.log(target,JSON.stringify({examined:r.examined,recipe:r.recipes[0]?.recipe,G:r.recipes[0]?.prediction.G_Pa,mc:r.recipes[0]?.mc.preparation90_Pa,p:r.recipes[0]?.mc.target_fraction_given_valid,outside:r.recipes[0]?.mc.outside_training_fraction}));
 }
 const any=await mc.robustInverse({target_Pa:15000,tolerance_percent:10,family:'any',ion:'Ca',supported_only:true},{n:5000,seed:20261014});fs.writeFileSync(out+'/search_any.json',JSON.stringify(any,null,2));console.log('any',any.examined,any.recipes.map(r=>r.prediction.family));
})();
