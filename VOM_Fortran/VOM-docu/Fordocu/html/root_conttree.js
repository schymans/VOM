  var tree;
  function treeInit()
  {
    tree = new YAHOO.widget.TreeView("treeDiv1");
    tree.setDynamicLoad(loadDataForNode);
    tree.subscribe("collapse", function(node){ 
       // -- remove the dynamic loaded child data after collapse of node --
       tree.removeChildren(node);
    }); 
    var root = tree.getRoot();
var myobj = { label: "coreprog.f90", id: "coreprog__f90", href: "", target:"basefrm" };
var tmpNode1 = new YAHOO.widget.TextNode(myobj, root, false);
var myobj = { label: "modules.f90", id: "modules__f90", href: "", target:"basefrm" };
var tmpNode3 = new YAHOO.widget.TextNode(myobj, root, false);
var myobj = { label: "readdata.f90", id: "readdata__f90", href: "", target:"basefrm" };
var tmpNode9 = new YAHOO.widget.TextNode(myobj, root, false);
var myobj = { label: "sample.f90", id: "sample__f90", href: "", target:"basefrm" };
var tmpNode11 = new YAHOO.widget.TextNode(myobj, root, false);
var myobj = { label: "sce.f90", id: "sce__f90", href: "", target:"basefrm" };
var tmpNode15 = new YAHOO.widget.TextNode(myobj, root, false);
var myobj = { label: "transpmodel.f90", id: "transpmodel__f90", href: "", target:"basefrm" };
var tmpNode33 = new YAHOO.widget.TextNode(myobj, root, false);
var myobj = { label: "watbal.f90", id: "watbal__f90", href: "", target:"basefrm" };
var tmpNode65 = new YAHOO.widget.TextNode(myobj, root, false);
tree.draw();
}

function loadDataForNode(node, onCompleteCallback)
{
   var id= node.data.id;
   // -- load the data, use dynamic functions defined in *_conttree.js --
    window[id](node,onCompleteCallback);
   // --scroll expanded node to top of view (if possible) 
   document.getElementById(node.contentElId).scrollIntoView(true);
}

