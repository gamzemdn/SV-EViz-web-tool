(function () {
    function setUploadStatus(filename) {
        const target = document.getElementById("tru-ref-upload-live-status");

        if (!target) {
            return;
        }

        target.innerHTML = `
            <div class="alert alert-warning" style="margin-top:6px; font-size:13px;">
                <strong>📤 Reference genome upload started...</strong>
                <div style="font-size:13px; margin-top:4px;">File: ${filename}</div>
                <div style="font-size:12px; margin-top:4px;">
                    Large FASTA files may take several minutes. Please do not close this tab.
                </div>
                <div class="progress" style="height:8px; margin-top:8px;">
                    <div class="progress-bar progress-bar-striped progress-bar-animated"
                         role="progressbar"
                         style="width:100%;">
                    </div>
                </div>
            </div>
        `;
    }

    function attachUploadListener() {
        const uploadDiv = document.getElementById("tru-ref-upload");

        if (!uploadDiv) {
            return;
        }

        const input = uploadDiv.querySelector('input[type="file"]');

        if (!input) {
            return;
        }

        if (input.dataset.truRefListenerAttached === "true") {
            return;
        }

        input.dataset.truRefListenerAttached = "true";

        input.addEventListener("change", function (event) {
            const file = event.target.files && event.target.files[0];

            if (file) {
                setUploadStatus(file.name);
            }
        });
    }

    setInterval(attachUploadListener, 1000);
})();
