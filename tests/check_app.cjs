// Offline application smoke checks and cross-language numerical fixture export.
const fs=require('fs'),path=require('path'),vm=require('vm'),assert=require('assert');
const R=path.resolve(__dirname,'..');
const html=fs.readFileSync(path.join(R,'index.html'),'utf8');
const code=html.match(/<script>([\s\S]*?)<\/script>/)[1];
const defaults={s0:'1.7',s1:'1.7',s2:'1.7',f0:'.65',f1:'.35',f2:'0',ion:'Ca',dose:'.06757',mcStock:'.1',mcVolume:'2',mcIon:'5',mcN:'100',mcSeed:'20261014',target:'15000',tolerance:'10',family:'ALG+TOCNF',searchIon:'Ca'};
const nodes={};const get=id=>nodes[id]??(nodes[id]={value:defaults[id]??'',checked:true,textContent:'',innerHTML:'',disabled:false,scrollIntoView(){}});
const context={console,setTimeout,Blob,URL,document:{getElementById:get,createElement:()=>({click(){}})}};
vm.createContext(context);vm.runInContext(code,context);
(async()=>{
 assert(get('prediction').textContent.includes('kPa'));
 const p=vm.runInContext('lastPrediction',context);assert(Math.abs(p.G_Pa-13650.0197)<.01);
 await get('simulate').onclick();assert(get('mcResult').innerHTML.includes('Conditional 5th'));
 assert(get('mcProgress').textContent.includes('100 valid'));
 get('mcmc').onclick();assert(get('mcmcResult').textContent.includes('G_Pa'));
 get('f0').value='0.5';get('predict').onclick();assert(!get('prediction').textContent.includes('kPa'));
 get('f0').value='0,65';get('predict').onclick();assert(get('prediction').textContent.includes('kPa'));
 await get('search').onclick();const search=vm.runInContext('lastSearch',context);assert(search.recipes.length>0);assert(search.examined===7182);
 vm.runInContext('useRecipe(0)',context);assert(get('prediction').textContent.includes('kPa'));
 const out={standalone_initialisation:true,forward_prediction:true,monte_carlo_handler:true,mcmc_handler:true,invalid_fraction_rejected:true,decimal_comma:true,inverse_handler:true,recipe_transfer:true,nominal_fixture_Pa:p.G_Pa,network_requests_required:false};
 fs.mkdirSync(path.join(R,'results'),{recursive:true});fs.writeFileSync(path.join(R,'results/app_verification.json'),JSON.stringify(out,null,2));console.log(JSON.stringify(out,null,2));
})().catch(e=>{console.error(e);process.exit(1)});
