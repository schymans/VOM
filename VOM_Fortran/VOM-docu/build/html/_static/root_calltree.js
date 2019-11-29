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
var myobj = { label: "vom", id: "vom", href: "vom.html", target:"basefrm" };
var tmpNode1 = new YAHOO.widget.TextNode(myobj, root, false);
tree.draw();
}

function loadDataForNode(node, onCompleteCallback)
{
   var id= node.data.id;
   // -- load the data, use dynamic functions defined in *_calltree.js --
    window[id](node,onCompleteCallback);
   // --scroll expanded node to top of view (if possible) 
   document.getElementById(node.contentElId).scrollIntoView(true);
}

