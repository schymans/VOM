function modules__f90(node, onCompleteCallback)
{
   var myobj = { label: "vom_file_mod", id: "vom_file_mod", href: "vom_file_mod.html", target:"basefrm" };
   var tmpNode4 = new YAHOO.widget.TextNode(myobj, node, false);
   tmpNode4.isLeaf = true; 
   var myobj = { label: "vegmod", id: "vegmod", href: "vegmod.html", target:"basefrm" };
   var tmpNode5 = new YAHOO.widget.TextNode(myobj, node, false);
   tmpNode5.isLeaf = true; 
   var myobj = { label: "watmod", id: "watmod", href: "watmod.html", target:"basefrm" };
   var tmpNode6 = new YAHOO.widget.TextNode(myobj, node, false);
   tmpNode6.isLeaf = true; 
   var myobj = { label: "vom_vegwat_mod", id: "vom_vegwat_mod", href: "vom_vegwat_mod.html", target:"basefrm" };
   var tmpNode7 = new YAHOO.widget.TextNode(myobj, node, false);
   tmpNode7.isLeaf = true; 
   var myobj = { label: "vom_sce_mod", id: "vom_sce_mod", href: "vom_sce_mod.html", target:"basefrm" };
   var tmpNode8 = new YAHOO.widget.TextNode(myobj, node, false);
   tmpNode8.isLeaf = true; 
     // notify the TreeView component when data load is complete
     onCompleteCallback();
}
