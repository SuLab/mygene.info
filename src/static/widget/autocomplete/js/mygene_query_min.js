mygene={url_root:"http://mygene.info/widget/autocomplete/",input_selector:"input.mygene_query_target",default_select_callback_name:"mygene_query_select_callback",loadfile:function(d,a,c){if("js"==a){var b=document.createElement("script");b.setAttribute("type","text/javascript");b.setAttribute("src",d)}else"css"==a&&(b=document.createElement("link"),b.setAttribute("rel","stylesheet"),b.setAttribute("type","text/css"),b.setAttribute("href",d));"undefined"!=typeof b&&(b.readyState?b.onreadystatechange=
function(){"complete"!=this.readyState&&"loaded"!=this.readyState||c()}:b.onload=c);(document.getElementsByTagName("head")[0]||document.getElementsByTagName("body")[0]).appendChild(b)},add_css:function(d){var a=document.createElement("style");a.type="text/css";a.styleSheet?a.styleSheet.cssText=d:a.appendChild(document.createTextNode(d));document.getElementsByTagName("head")[0].appendChild(a)},genequery:function(d){var a=window.$||window.JQuery,c=a(mygene.input_selector);c.autocomplete({source:function(b,
c){a.ajax({url:"http://mygene.info/v2/query",dataType:"jsonp",jsonp:"callback",data:{q:b.term,species:"human",size:20},success:function(b){0<b.total?c(a.map(b.hits,function(b){return a.extend(b,{label:b.symbol+": "+b.name,id:b._id,value:b.symbol})})):c([{label:"no matched gene found.",value:""}])}})},minLength:2,select:d||mygene.default_select_callback,open:function(){a(this).removeClass("ui-corner-all").addClass("ui-corner-top")},close:function(){a(this).removeClass("ui-corner-top").addClass("ui-corner-all")}});
void 0===c.attr("title")&&c.attr("title","Powered by mygene.info")},default_select_callback:function(d,a){alert(a.item?"Selected: "+a.item.label+"("+a.item.id+")":"Nothing selected, input was "+this.value)}};
mygene_init=function(){function d(){var c=window.$||window.JQuery;void 0===c.ui||"1.8.21"!==c.ui.version?mygene.loadfile(mygene.url_root+"js/jquery-ui-1.8.21.custom.min.js","js",a):a()}function a(){mygene.genequery(window[mygene.default_select_callback_name]);mygene.loadfile(mygene.url_root+"css/ui-lightness/jquery-ui-1.8.21.custom.css","css");mygene.add_css('.ui-autocomplete-loading { background: white url("'+mygene.url_root+'img/ui-anim_basic_16x16.gif") right center no-repeat; }')}(function(){var a=
window.$||window.JQuery;void 0===a||"1.7.2"!==a.fn.jquery?mygene.loadfile(mygene.url_root+"js/jquery-1.7.2.min.js","js",d):d()})()};mygene_init();