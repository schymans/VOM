function sce(node, onCompleteCallback)
{
   var myobj = { label: "sce_init", id: "sce_init", href: "sce_init.html", target:"basefrm" };
   var tmpNode33 = new YAHOO.widget.TextNode(myobj, node, false);
   var myobj = { label: "initialseed", id: "initialseed", href: "initialseed.html", target:"basefrm" };
   var tmpNode34 = new YAHOO.widget.TextNode(myobj, node, false);
   var myobj = { label: "sortcomp", id: "sortcomp", href: "sortcomp.html", target:"basefrm" };
   var tmpNode68 = new YAHOO.widget.TextNode(myobj, node, false);
   var myobj = { label: "writepars", id: "writepars", href: "writepars.html", target:"basefrm" };
   var tmpNode70 = new YAHOO.widget.TextNode(myobj, node, false);
   var myobj = { label: "run_cce", id: "run_cce", href: "run_cce.html", target:"basefrm" };
   var tmpNode71 = new YAHOO.widget.TextNode(myobj, node, false);
   var myobj = { label: "ck_success", id: "ck_success", href: "ck_success.html", target:"basefrm" };
   var tmpNode108 = new YAHOO.widget.TextNode(myobj, node, false);
     // notify the TreeView component when data load is complete
     onCompleteCallback();
}
