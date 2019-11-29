function run_cce(node, onCompleteCallback)
{
   var myobj = { label: "omp_set_num_threads", id: "omp_set_num_threads", href: "", target:"basefrm" };
   var tmpNode58 = new YAHOO.widget.TextNode(myobj, node, false);
   tmpNode58.isLeaf = true; 
   var myobj = { label: "vom_dealloc", id: "vom_dealloc", href: "vom_dealloc.html", target:"basefrm" };
   var tmpNode59 = new YAHOO.widget.TextNode(myobj, node, false);
   tmpNode59.isLeaf = true; 
   var myobj = { label: "cce", id: "cce", href: "cce.html", target:"basefrm" };
   var tmpNode60 = new YAHOO.widget.TextNode(myobj, node, false);
     // notify the TreeView component when data load is complete
     onCompleteCallback();
}
