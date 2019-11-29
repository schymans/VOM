function vom(node, onCompleteCallback)
{
   var myobj = { label: "read_commandline", id: "read_commandline", href: "read_commandline.html", target:"basefrm" };
   var tmpNode2 = new YAHOO.widget.TextNode(myobj, node, false);
   var myobj = { label: "read_shufflepar", id: "read_shufflepar", href: "read_shufflepar.html", target:"basefrm" };
   var tmpNode4 = new YAHOO.widget.TextNode(myobj, node, false);
   tmpNode4.isLeaf = true; 
   var myobj = { label: "transpmodel_init_once", id: "transpmodel_init_once", href: "transpmodel_init_once.html", target:"basefrm" };
   var tmpNode5 = new YAHOO.widget.TextNode(myobj, node, false);
   var myobj = { label: "sce_main", id: "sce_main", href: "sce_main.html", target:"basefrm" };
   var tmpNode15 = new YAHOO.widget.TextNode(myobj, node, false);
   var myobj = { label: "transpmodel", id: "transpmodel", href: "transpmodel.html", target:"basefrm" };
   var tmpNode135 = new YAHOO.widget.TextNode(myobj, node, false);
   var myobj = { label: "random_samples", id: "random_samples", href: "random_samples.html", target:"basefrm" };
   var tmpNode165 = new YAHOO.widget.TextNode(myobj, node, false);
     // notify the TreeView component when data load is complete
     onCompleteCallback();
}
